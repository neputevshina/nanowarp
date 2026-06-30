package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
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

	fft         *fourier.FFT
	arm         []bool  // PGHI done mask
	norm, wgain float64 // Global normalization factor and grain-only normalization factor
	heap        hp      // PGHI heap

	a wbufs
}

type wbufs struct {
	S, Mid, M, P  []float64      // Scratch buffers
	Ph            []float64      `size:"nbins"` // Current phase
	Past          []float64      // Phase accumulators
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
		lah:   olap + 1,
	}
	a := &n.a

	makeslices(a, n.nbins, nfft, nch, n.lah)
	n.arm = make([]bool, n.nbins)

	s := func(w []float64) []float64 {
		return w[:nbuf]
	}
	blackmanHarris(s(a.W))
	n.heap = make(hp, n.nbins)

	// FIXME Destination is the first argument by convention.
	windowDx(s(a.W), s(a.Wd))
	windowT(s(a.W), s(a.Wt))

	copy(s(a.Wr), s(a.W))
	slices.Reverse(s(a.Wr))
	n.wgain = windowGain(n.a.W)
	n.norm = float64(nfft) * float64(n.olap) * n.osamp * n.wgain

	n.fft = fourier.NewFFT(nfft)
	n.heap = make(hp, n.lah*n.nbins) // 2 for future and past.

	return
}

func (n *warper) process3(in [][]float64, out [][]float64, coeffs, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	get := func() [][]float64 { return make2[float64](len(in), n.nfft) }
	grainbuf := get()
	ingrain := get()

	nch := len(in)

	lastone := 0
	for j := -n.nbuf; ; j += n.hop {
		if j > len(out[0])-n.nbuf {
			break
		}
		fi := func(j int) int { return int(phasor[max(0, j)]) }
		i := fi(j)
		c := 1 / coeffs[max(0, j)]

		for ch := range nch {
			copy(ingrain[ch][max(0, n.nbuf/2-i):],
				in[ch][max(0, (i-n.nbuf/2)):clamp(0, len(in[ch]), i+n.nbuf/2)])
		}

		ph, df, _ := n.advance(ingrain, abs(c), c == 1)
		n.synthesize(grainbuf, ph, df)

		// Cut pre-echo in transient regions.
		d := j - lastone
		if c == 1 {
			lastone = j
		} else if d < n.nbuf/2 {
			z := func(grain []float64) {
				rr := grain[max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grain[:n.nbuf/2-d-n.hop], 0)
			}
			for ch := range nch {
				z(grainbuf[ch])
			}
		}

		if n.root.opts.Onsets && c != 1 {
			for ch := range nch {
				clear(grainbuf[ch])
			}
		}

		for ch := range nch {
			g := out[ch][max(0, j-n.nbuf/2):clamp(0, len(out[0]), j+n.nbuf/2)]
			add(g, grainbuf[ch][clamp(0, n.nbuf, -j):])
		}
	}
}

// advance constructs the next frame of the output.
func (n *warper) advance(ingrain [][]float64, stretch float64, reset bool) (ph []complex128, c [][]complex128, m []float64) {
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

	// Encode stereo phase differences and stretch mid only, keep original magnitudes.
	// NB: Phase difference in polar coordinates is complex division in cartesian.
	//     Phase sum is conversely a multiply.
	//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
	//
	// See Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
	// real-time realization for a commercial music product
	for w := range a.X {
		a.M[w] = mag(a.X[w])
		a.Y[w] = safediv(a.X[w], complex(a.M[w], 0))
		for ch := range len(ingrain) {
			// BREAKING: significant null test fail with commit 9e8e1747 due to numerical
			// instablility of a/b*b.
			//
			// Previously, a.C was bypassed on resets, now it is divided and then
			// (in (*warper).synthesize) multiplied back.
			// Required for GLA, which in other case would intoduce insignificant
			// null test fail (like one introduced in commit c7d6ddbb).
			//
			// Transients are still more or less left intact.
			a.C[ch][w] /= a.Y[w]
		}
	}

	// Bypass and refill the past on a phase reset.
	if reset {
		for w := range a.Y {
			a.Ph[w] = cmplx.Phase(a.X[w])
		}
		goto skip
	}

	// Calculate derivative windows.
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// Calculate partial derivatives of the phase.
	//
	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms
	//
	// TODO Probably it will be more numerically stable to limit the phase accuum to
	// 0..1 and scale back to -π..π range at the poltocar conversion.
	// fadv must return 0..1 accordingly, simply defer π multiplication until the end.
	for w := range a.X {
		a.Fadv[w] = princarg(fadv(a.X[:n.nbins], a.Xt, stretch, w))
		a.Tadv[w] = tadv(a.X[:n.nbins], a.Xd, float64(n.nfft/n.hop), w)
	}

	// Prepare heap (priority queue).
	n.heap = n.heap[:n.nbins]
	fill(n.arm, true)
	for w := range a.X {
		n.heap[w] = heaptriple{a.P[w], w, -1}
	}
	heapInit(&n.heap)

	// Perform phase gradient heap integration.
	//
	// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
	// (https://arxiv.org/pdf/2202.07382)
	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[w] {
				a.Ph[w] = princarg(a.Past[w]) + a.Tadv[w]
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, 0})
			}
		case 0:
			if w >= 1 && n.arm[w-1] {
				a.Ph[w-1] = a.Ph[w] - a.Fadv[w-1]
				n.arm[w-1] = false
				heapPush(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				a.Ph[w+1] = a.Ph[w] + a.Fadv[w+1]
				n.arm[w+1] = false
				heapPush(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
			}
		}
	}

	// Add stereo phase differences back through multiplication
	// and update past phases.
	for w := range a.Ph {
		a.Y[w] = cmplx.Rect(1, a.Ph[w])
	}
skip:
	copy(a.P, a.M)
	copy(a.Past, a.Ph)

	return a.Y, a.C, a.M
}

func (n *warper) synthesize(outgrain [][]float64, ph []complex128, c [][]complex128) {
	a := &n.a
	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.Wr)
		copy(out, a.S)
	}
	for ch := range outgrain {
		mul(c[ch], ph)
		defft(outgrain[ch], c[ch])
	}
}

// fadv calculates the partial derivative of the phase with respect
// to frequency using time-frequency reassignment.
func fadv(x, xt []complex128, stretch float64, w int) float64 {
	if mag(x[w]) == 0 {
		return 0
	}
	// NOTE Try len(x)-1 instead. Sounds worse on my $4 speakers.
	// FIXME Works ONLY with nbuf=4096, nfft=8192 (oversampling 2).
	return -real(xt[w]/x[w])/float64(len(x))*math.Pi*stretch - math.Pi/2
}

// tadv calculates the partial derivative of the phase with respect
// to time using time-frequency reassignment.
func tadv(x, xd []complex128, olap float64, w int) float64 {
	if mag(x[w]) < 1e-6 {
		return 0
	}
	return (math.Pi*float64(w) + imag(xd[w]/x[w])) / (olap / 2)
}
