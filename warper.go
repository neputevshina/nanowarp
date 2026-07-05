package nanowarp

import (
	"fmt"
	"math"
	"math/bits"
	"math/cmplx"
	"os"
	"slices"

	"gonum.org/v1/gonum/cmplxs"
	"gonum.org/v1/gonum/dsp/fourier"
	"gonum.org/v1/gonum/floats"
)

type warper struct {
	nfft  int     // DFT size, a power of 2
	nbuf  int     // Effective window size, nbuf<nfft
	hop   int     // Window output hop size
	lah   int     // Non-causal PGHI lookahead in frames (hop sizes)
	nbins int     // nfft/2+1, Number of DFT bins
	olap  int     // nbuf/hop, Window ovelap
	osamp float64 // nfft/nbuf, Zero-padding ratio

	root *Nanowarp

	fft *fourier.FFT

	// PGHI-related
	arm    []bool    // Done mask
	heap   hp        // Heap
	arrows [][2]int  // Integration directions
	trace  []float64 // Ridge accumulator
	ridges []uint    // Extracted ridges

	norm, wgain float64 // Global normalization factor and grain-only normalization factor

	a wbufs
}

type wbufs struct {
	S, Mid        []float64      // Scratch buffers
	M, P          []float64      `size:"nbins"` // Current and previous magnitude
	Ph            []float64      `size:"nbins"` // Current phase
	Past          []float64      // Phase accumulator
	Fadv, Tadv    []float64      // Partial derivatives
	W, Wr, Wd, Wt []float64      // Window functions
	X, Y, Xd, Xt  []complex128   // Complex spectra
	C             [][]complex128 // Channel differences
}

func warperNew(nbuf, osamp, olap, nch int, nanowarp *Nanowarp) (n *warper) {
	// FIXME Only 2x oversampling works, no more, no less.
	nfft := nextpow2(nbuf * osamp)
	n = &warper{
		nfft:  nfft,
		nbins: nfft/2 + 1,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		olap:  olap,
		osamp: float64(osamp),
		root:  nanowarp,
		lah:   3 * olap,
	}
	a := &n.a

	makeslices(a, n.nbins, nfft, nch, n.lah)
	n.arm = make([]bool, n.nbins)

	s := func(w []float64) []float64 {
		return w[:nbuf]
	}
	blackmanHarris(s(a.W))
	// hann(s(a.W))
	n.heap = make(hp, n.nbins)

	windowDx(s(a.Wd), s(a.W))
	windowT(s(a.Wt), s(a.W))

	copy(s(a.Wr), s(a.W))
	slices.Reverse(s(a.Wr))
	n.wgain = windowGain(n.a.W)
	n.norm = float64(nfft) * float64(n.olap) * n.osamp * n.wgain

	n.fft = fourier.NewFFT(nfft)
	n.heap = make(hp, 2*n.nbins) // 2 for future and past.
	n.arrows = make([][2]int, 0, n.nbins)
	n.ridges = make([]uint, n.nbins)
	n.trace = make([]float64, n.nbins)

	return
}

func (n *warper) process6(in [][]float64, out [][]float64, coeffs, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	get := func() [][]float64 { return make2[float64](len(in), n.nfft) }
	nch := len(in)

	lead := get()
	grain := get()
	crop := make([][]float64, nch)

	lastone := 0
	for j := -n.nbuf / 2; j < len(out[0])-1+n.nbuf/2; j += n.hop {
		bounds := func(i int) int { return clamp(0, len(out[0])-1, i) }
		fi := func(j int) int { return int(phasor[bounds(j)]) }
		i := fi(j)
		c := 1 / coeffs[bounds(j)]

		for ch := range nch {
			cr := in[ch][max(0, (i-n.nbuf/2)):clamp(0, len(in[ch]), i+n.nbuf/2)]
			if i < n.nbuf/2 {
				copy(lead[ch][max(0, n.nbuf/2-i):], cr)
				crop[ch] = lead[ch]
			} else {
				crop[ch] = cr
			}
		}

		normal, diff, _ := n.advance(crop, abs(c), n.root.opts.Quality >= 0 && c == 1)
		n.synthesize(grain, normal, diff)

		d := j - lastone
		if c == 1 {
			lastone = j
		}
		for ch := range nch {
			// Cut pre-echo in transient regions.
			if c != 1 && d < n.nbuf/2 {
				rr := grain[ch][max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grain[ch][:n.nbuf/2-d-n.hop], 0)
			}

			if n.root.opts.Onsets && c != 1 {
				clear(grain[ch])
			}

			g := out[ch][max(0, j-n.nbuf/2):bounds(j+n.nbuf/2)]
			add(g, grain[ch][clamp(0, n.nbuf, -j):])
		}
	}
}

// advance constructs the next frame of the output.
func (n *warper) advance(ingrain [][]float64, stretch float64, reset bool) (normal []complex128, diff [][]complex128, mag []float64) {
	a := &n.a
	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	clear(a.Mid)
	for ch := range ingrain {
		add(a.Mid, ingrain[ch])
		enfft(a.C[ch], a.W, ingrain[ch])
	}
	enfft(a.X, a.W, a.Mid)

	cmplxs.Abs(a.M, a.X)
	arr := n.pghiarrows(a.P, a.M, n.arm, n.arrows, n.ridges)

	// Calculate time of life for each ridge.
	//
	// Note that trajectories are calculated with a 1-frame lag.
	// Still good for our purposes.
	trace := n.trace
	for w, v := range n.ridges {
		p := boolfloat(bits.OnesCount(v&0b111000) >= 2)
		trace[w] = trace[w]*p + p
	}
	l := -1
	for i := range trace {
		if l > 0 && trace[i] == 0 {
			fill(trace[l:i], slices.Max(trace[l:i]))
		}
		if trace[i] != 0 {
			l = i
		}
	}
	if l > 0 {
		fill(trace[l:], slices.Max(trace[l:]))
	}

	// Bypass on a phase reset.
	if reset {
		for w := range a.Y {
			a.Ph[w] = cmplx.Phase(a.X[w])
		}
		fill(a.Y, 1)
		goto skip
	}

	// Encode stereo phase differences and stretch mid only, keep original magnitudes.
	// NB: Phase difference in polar coordinates is complex division in cartesian.
	//     Phase sum is conversely a multiply.
	//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
	//
	// See Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
	// real-time realization for a commercial music product
	for w := range a.X {
		a.Y[w] = safediv(a.X[w], complex(a.M[w], 0))
	}
	for ch := range len(ingrain) {
		cmplxs.Div(a.C[ch], a.Y)
	}

	// Calculate partial derivatives of the phase.
	//
	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms
	//
	// TODO Probably it will be more numerically stable to limit the phase accuum to
	// 0..1 and scale back to -π..π range at the poltocar conversion.
	// fadv must return 0..1 accordingly, simply defer π multiplication until the end.
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)
	for w := range a.X {
		a.Fadv[w] = princarg(fadv(a.X[:n.nbins], a.Xt, stretch, w))
		a.Tadv[w] = tadv(a.X[:n.nbins], a.Xd, float64(n.nfft/n.hop), w)
	}

	n.pghiintegrate(arr, a.Fadv, a.Tadv, a.Ph, a.Past)

	// Add stereo phase differences back through multiplication.
	for w := range a.Ph {
		a.Y[w] = cmplx.Rect(1, a.Ph[w])
	}

skip:
	copy(a.P, a.M)
	copy(a.Past, a.Ph)

	return a.Y, a.C, a.M
}

func (n *warper) synthesize(outgrain [][]float64, normal []complex128, diff [][]complex128) {
	a := &n.a
	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		floats.Scale(1/n.norm, a.S)
		mul(a.S, a.Wr)
		copy(out, a.S)
	}

	for ch := range outgrain {
		mul(diff[ch], normal)
		defft(outgrain[ch], diff[ch])
	}
}

// pghiarrows calculates integration directions from a pair of magnitude
// buffers using priority queue.
//
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) pghiarrows(P, M []float64, arm []bool, arrows [][2]int, ridges []uint) [][2]int {
	n.heap = n.heap[:n.nbins]
	fill(n.arm, true)
	for w := range P {
		n.heap[w] = heaptriple{P[w], w, -1}
		// Ridges are encoded as 3-bit mask: 1 right, 2 down, 4 up
		ridges[w] <<= 3
	}
	heapInit(&n.heap)
	arrows = arrows[:0]

	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		w := h.w
		switch h.t {
		case -1:
			if arm[w] {
				arrows = append(arrows, [2]int{w, 'p'})
				ridges[w] |= 1 << 3
				arm[w] = false
				heapPush(&n.heap, heaptriple{M[w], w, 0})
			}
		case 0:
			if w >= 1 && arm[w-1] {
				arrows = append(arrows, [2]int{w, 'd'})
				ridges[w-1] |= 2
				arm[w-1] = false
				heapPush(&n.heap, heaptriple{M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && arm[w+1] {
				arrows = append(arrows, [2]int{w, 'u'})
				ridges[w+1] |= 4
				arm[w+1] = false
				heapPush(&n.heap, heaptriple{M[w+1], w + 1, 0})
			}
		}
	}
	return arrows
}

// pghiintegrate integrates partial derivatives of phase using directions
// obtained from pghiarrows.
func (n *warper) pghiintegrate(arrows [][2]int, Fadv, Tadv, Ph, Past []float64) {
	for _, e := range arrows {
		switch e[1] {
		case 'p':
			Ph[e[0]] = princarg(Past[e[0]]) + Tadv[e[0]]
		case 'd':
			Ph[e[0]-1] = Ph[e[0]] - Fadv[e[0]-1]
		case 'u':
			Ph[e[0]+1] = Ph[e[0]] + Fadv[e[0]+1]
		}
	}
}

// fadv calculates the partial derivative of the phase with respect
// to frequency using time-frequency reassignment.
func fadv(x, xt []complex128, stretch float64, w int) float64 {
	if cmplx.Abs(x[w]) == 0 {
		return 0
	}
	// NOTE Try len(x)-1 instead. Sounds worse on my $4 speakers.
	return -real(xt[w]/x[w])/float64(len(x))*math.Pi*stretch - math.Pi/2
}

// tadv calculates the partial derivative of the phase with respect
// to time using time-frequency reassignment.
func tadv(x, xd []complex128, olap float64, w int) float64 {
	if cmplx.Abs(x[w]) < 1e-6 {
		return 0
	}
	return (math.Pi*float64(w) + imag(xd[w]/x[w])) / (olap / 2)
}
