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
	nfft  int
	nbuf  int
	nbins int
	hop   int
	olap  int
	root  *Nanowarp

	fft  *fourier.FFT
	arm  []bool
	norm float64
	heap hp
	img  [][]float64

	a wbufs
}
type wbufs struct {
	S, Mid, M, P, Phase, Pphase        []float64 // Scratch buffers
	Fadv, Pfadv, Ldiff, Rdiff          []float64
	W, Wr, Wd, Wt                      []float64    // Window functions
	X, Y, Xd, Xt, L, Ld, R, Rd, Lo, Ro []complex128 // Complex spectra
}

func warperNew(nbuf int, nanowarp *Nanowarp) (n *warper) {
	nfft := nextpow2(nbuf) * 2
	nbins := nfft/2 + 1
	olap := 4
	if nanowarp.opts.Quality >= 1 {
		olap = 8
	}
	n = &warper{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		olap:  olap,
		root:  nanowarp,
	}
	a := &n.a

	makeslices(a, nbins, nfft)

	// Exceptions
	a.Phase = make([]float64, nbins)
	n.arm = make([]bool, nbins)

	blackmanHarris(a.W[:nbuf])

	windowDx(a.W[:nbuf], a.Wd[:nbuf])
	windowT(a.W[:nbuf], a.Wt[:nbuf])
	copy(a.Wr[:nbuf], a.W[:nbuf])
	slices.Reverse(a.Wr[:nbuf])
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * windowGain(n.a.W)

	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *warper) process3(lin, rin, lout, rout []float64, coeffs, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	get := func() []float64 { return make([]float64, n.nfft) }
	lgrainbuf, rgrainbuf := get(), get()
	lingrain, ringrain := get(), get()

	lastone := 0
	for j := -n.nbuf; ; j += n.hop {
		if j > len(lout)-n.nbuf {
			break
		}
		i := int(phasor[max(0, j)] - float64(n.nbuf/2))

		if i > len(lin)-n.nbuf {
			break
		}
		clear(lingrain)
		clear(ringrain)
		copy(lingrain[max(0, -i):], lin[max(0, i):clamp(0, len(lin), i+n.nbuf)])
		copy(ringrain[max(0, -i):], rin[max(0, i):clamp(0, len(lin), i+n.nbuf)])
		c := 1.
		if j > 0 {
			c = 1 / coeffs[j]
			if c != c || math.IsInf(c, 0) {
				c = 1
			}
		}
		tensec := 10 * n.root.fs
		if n.root.opts.Quality == -1 && j%tensec < (j-n.hop)%tensec {
			c = 1
		}

		n.advance(lingrain, ringrain, lgrainbuf, rgrainbuf, abs(c), c == 1)

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
			z(lgrainbuf)
			z(rgrainbuf)
		}

		if n.root.opts.Onsets && c != 1 {
			clear(lgrainbuf)
			clear(rgrainbuf)
		}

		loutgrain := lout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		routgrain := rout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]

		add(loutgrain, lgrainbuf[clamp(0, n.nbuf, -j):])
		add(routgrain, rgrainbuf[clamp(0, n.nbuf, -j):])
	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(lingrain, ringrain, loutgrain, routgrain []float64, stretch float64, reset bool) {
	a := &n.a
	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	for i := range lingrain {
		a.Mid[i] = lingrain[i] + ringrain[i]
	}

	enfft(a.L, a.W, lingrain)
	enfft(a.R, a.W, ringrain)

	enfft(a.X, a.W, a.Mid)
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	// TODO Works ONLY with 2x FFT oversampling.
	olap := float64(n.nbuf / n.hop)
	osampc := 2.
	tadv := gettadv(a.X[:n.nbins], a.Xd, osampc, olap)
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch)

	// Force phase reset.
	// TODO Move out to process3, PV doesn't have to handle phase resets.
	if reset {
		for w := range a.Phase {
			a.Pphase[w] = cmplx.Phase(a.X[w])
		}
		copy(a.P, a.M)
		copy(a.Lo, a.L)
		copy(a.Ro, a.R)
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

		a.L[w] /= p
		a.R[w] /= p
	}

	n.heap = make(hp, n.nbins)
	clear(n.arm)

	for j := range a.X {
		n.arm[j] = true
		n.heap[j] = heaptriple{a.P[j], j, -1}
	}
	heapInit(&n.heap)

	// Do PGHI.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap).(heaptriple)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[w] {
				adv := tadv(w)
				a.Phase[w] = a.Pphase[w] + adv
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, 0})
			}
		case 0:
			if w > 1 && n.arm[w-1] {
				adv := -a.Fadv[w-1]
				a.Phase[w-1] = a.Phase[w] + adv
				n.arm[w-1] = false
				heapPush(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				adv := a.Fadv[w+1]
				a.Phase[w+1] = a.Phase[w] + adv
				n.arm[w+1] = false
				heapPush(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
			}
		}
	}

	copy(a.P, a.M)
	for w := range a.Phase {
		// Add stereo phase differences back through multiplication.
		a.Lo[w] = a.L[w] * cmplx.Rect(1, a.Phase[w])
		a.Ro[w] = a.R[w] * cmplx.Rect(1, a.Phase[w])
		a.Pphase[w] = princarg(a.Phase[w])
	}
	goto skip
skip:
	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.Wr)
		copy(out, a.S)
	}
	defft(loutgrain, a.Lo)
	defft(routgrain, a.Ro)
}
