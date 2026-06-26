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
	arm         [][]bool // PGHI done mask
	norm, wgain float64  // Global normalization factor and grain-only normalization factor
	heap        hp       // PGHI heap

	a wbufs
}

type wbufs struct {
	S, Mid, M, P, F    []float64      // Scratch buffers
	Ph                 []float64      `size:"nbins"` // Current phase
	Past               []float64      // Phase accumulators
	Fadv, Tadv         []float64      // Partial derivatives
	W, Wr, Wd, Wt, Wdt []float64      // Window functions
	X, Y, Xd, Xt       []complex128   // Complex spectra
	C, Co              [][]complex128 // Channels
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

	n.arm = make([][]bool, n.lah)
	for i := range n.arm {
		n.arm[i] = make([]bool, n.nbins)
	}

	s := func(w []float64) []float64 {
		return w[:nbuf]
	}
	blackmanHarris(s(a.W))

	// FIXME Destination is the first argument by convention.
	windowDx(s(a.W), s(a.Wd))
	windowT(s(a.W), s(a.Wt))
	windowT(s(a.Wd), s(a.Wdt))

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

		cut := func(lingrain, lin []float64, i int) {
			copy(lingrain[max(0, n.nbuf/2-i):], lin[max(0, (i-n.nbuf/2)):clamp(0, len(lin), i+n.nbuf/2)])
		}
		for ch := range nch {
			cut(ingrain[ch], in[ch], i)
		}

		n.advance(ingrain, grainbuf, abs(c), c == 1)

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

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(ingrain, outgrain [][]float64, stretch float64, reset bool) {
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
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	// TODO Works ONLY with 2x FFT oversampling.
	olap := float64(n.nbuf / n.hop)
	osampc := 2.
	tadv := gettadv(a.X[:n.nbins], a.Xd, osampc*olap)
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch)

	if reset {
		for w := range a.Ph {
			a.Ph[w] = cmplx.Phase(a.X[w])
		}
		copy(a.Past, a.Ph)
		for ch := range len(ingrain) {
			copy(a.Co[ch], a.C[ch])
		}
		goto skip
	}

	for w := range a.X {
		// TODO Probably it will be more numerically stable to limit the phase accuum to
		// 0..1 and scale back to -π..π range at the poltocar conversion.
		// fadv must return 0..1 accordingly, simply defer π multiplication till the end.
		a.Fadv[w] = princarg(fadv(w))
	}

	for w := range a.X {
		// Encode stereo phase differences and stretch mid only, keep original magnitudes.
		// NB: Phase difference in polar coordinates is complex division in cartesian.
		//     Phase sum is conversely a multiply.
		//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
		//
		// See “Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
		// real-time realization for a commercial music product.”
		m := mag(a.X[w])
		p := a.X[w] / complex(m, 0)
		if m < 1e-6 {
			p = complex(1, 0)
		}
		a.M[w] = m
		for ch := range len(ingrain) {
			a.C[ch][w] /= p
		}
	}

	n.heap = make(hp, n.nbins)

	clear(n.arm[0])
	for w := range a.X {
		n.arm[0][w] = true
		n.heap[w] = heaptriple{a.P[w], w, -1}
	}

	heapInit(&n.heap)

	// Do PGHI.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[0][w] {
				adv := tadv(w)
				a.Ph[w] = a.Past[w] + adv
				n.arm[0][w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, 0})
			}
		case 0:
			if w > 1 && n.arm[0][w-1] {
				adv := -a.Fadv[w-1]
				a.Ph[w-1] = a.Ph[w] + adv
				n.arm[0][w-1] = false
				heapPush(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[0][w+1] {
				adv := a.Fadv[w+1]
				a.Ph[w+1] = a.Ph[w] + adv
				n.arm[0][w+1] = false
				heapPush(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
			}
		}
	}

	for w := range a.Ph {
		// Add stereo phase differences back through multiplication.
		phasor := cmplx.Rect(1, a.Ph[w])
		for ch := range ingrain {
			a.Co[ch][w] = a.C[ch][w] * phasor
		}
		a.Past[w] = princarg(a.Ph[w])
	}
	goto skip
skip:
	copy(a.P, a.M)
	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.Wr)
		copy(out, a.S)
	}
	for ch := range ingrain {
		defft(outgrain[ch], a.Co[ch])
	}
}

func getfadv(x, xt []complex128, stretch float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) == 0 {
			return 0
		}
		// NOTE Try len(x)-1 instead. Sounds worse on my $4 speakers.
		// FIXME Works ONLY with nbuf=4096, nfft=8192 (oversampling 2).
		return -real(xt[j]/x[j])/float64(len(x))*math.Pi*stretch - math.Pi/2
	}
}

func gettadv(x, xd []complex128, olap float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) < 1e-6 {
			return 0
		}
		return (math.Pi*float64(j) + imag(xd[j]/x[j])) / (olap / 2)
	}
}
