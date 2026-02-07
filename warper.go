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
	S, Mid, M, P, Phase, Pphase     []float64 // Scratch buffers
	Fadv, Pfadv, Ldiff, Rdiff       []float64
	W, Wr, Wd, Wt                   []float64    // Window functions
	X, Xd, Xt, L, Ld, R, Rd, Lo, Ro []complex128 // Complex spectra
	Pha, Px                         []complex128
}

func warperNew(nbuf int, nanowarp *Nanowarp) (n *warper) {
	nfft := nextpow2(nbuf) * 2
	nbins := nfft/2 + 1
	olap := 4
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

	// Exceptions.
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

func (n *warper) process3(lin, rin, lout, rout []float64, shift []float64, delay float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	get := func() []float64 { return make([]float64, n.nfft) }
	lgrainbuf := get()
	rgrainbuf := get()
	lingrain := get()
	ringrain := get()

	bayer := []float64{0.2, 0.6, 0.4, 0.8}
	for j := -n.nbuf; ; j += n.hop {
		di, df := math.Modf(delay)
		if j > len(lout)-n.nbuf {
			break
		}
		i := int(shift[clamp(0, len(lout)-1, j)])
		if i > len(lin)-n.nbuf {
			break
		}
		clear(lingrain)
		clear(ringrain)
		copy(lingrain[max(0, -i):], lin[max(0, i):min(len(lin), i+n.nbuf)])
		copy(ringrain[max(0, -i):], rin[max(0, i):min(len(lin), i+n.nbuf)])
		coeff := 1.
		if j > 0 {
			coeff = 1 / (shift[j] - shift[j-1])
			if coeff != coeff || math.IsInf(coeff, 0) {
				coeff = 1
			}
		}
		twosec := 2 * n.root.fs
		if n.root.opts.Raw && j%twosec < (j-n.hop)%twosec {
			coeff = 1
		}

		n.advance(lingrain, ringrain, lgrainbuf, rgrainbuf, abs(coeff))
		if n.root.opts.Onsets && coeff != 1 {
			clear(lgrainbuf)
			clear(rgrainbuf)
		}

		tj := j + int(di) + boolint(df > bayer[j%len(bayer)])
		loutgrain := lout[max(0, tj):clamp(0, len(lout), tj+n.nbuf)]
		routgrain := rout[max(0, tj):clamp(0, len(lout), tj+n.nbuf)]
		add(loutgrain, lgrainbuf[clamp(0, n.nbuf, -tj):])
		add(routgrain, rgrainbuf[clamp(0, n.nbuf, -tj):])
	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(lingrain, ringrain, loutgrain, routgrain []float64, stretch float64) {
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

	if n.root.opts.Diffadv {
		enfft(a.Ld, a.Wd, lingrain)
		enfft(a.Rd, a.Wd, ringrain)
	}

	enfft(a.X, a.W, a.Mid)
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	// TODO Works ONLY with 2x FFT oversampling.
	olap := float64(n.nbuf / n.hop)
	osampc := 2.
	tadv := gettadv(a.X[:n.nbins], a.Xd, osampc, olap)
	ltadv := gettadv(a.L[:n.nbins], a.Ld, osampc, olap)
	rtadv := gettadv(a.R[:n.nbins], a.Rd, osampc, olap)
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch)

	// Force reset on a no-op stretch, including detected transients.
	if stretch == 1 {
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
		if n.root.opts.Diffadv {
			a.Ldiff[w] = tadv(w) - ltadv(w)
			a.Rdiff[w] = tadv(w) - rtadv(w)
		} else {
			a.L[w] /= p
			a.R[w] /= p
		}
	}

	n.heap = make(hp, n.nbins)
	clear(n.arm)

	for j := range a.X {
		n.arm[j] = true
		n.heap[j] = heaptriple{a.P[j], j, -1}
	}
	heapInit(&n.heap)

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
		if n.root.opts.Diffadv {
			a.Lo[w] = cmplx.Rect(mag(a.L[w]), a.Phase[w]+a.Ldiff[w])
			a.Ro[w] = cmplx.Rect(mag(a.R[w]), a.Phase[w]+a.Rdiff[w])
		} else {
			// Add stereo phase differences back through multiplication.
			a.Lo[w] = a.L[w] * cmplx.Rect(1, a.Phase[w])
			a.Ro[w] = a.R[w] * cmplx.Rect(1, a.Phase[w])
		}
		a.Pphase[w] = princarg(a.Phase[w])
	}
	goto skip
skip:
	copy(a.Px, a.X)
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
