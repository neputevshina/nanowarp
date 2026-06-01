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

	fft         *fourier.FFT
	arm         []bool
	norm, wgain float64
	heap        hp
	img         [][]float64

	a wbufs
}
type wbufs struct {
	S, Mid, M, P               []float64 // Scratch buffers
	Phase, Rphase              []float64 // Current phases (advances)
	Pphase, Nphase             []float64 // Phase accumulators
	Fadv                       []float64
	W, Wr, Wd, Wt              []float64      // Window functions
	X, E, Xd, Xt, L, R, Lo, Ro []complex128   // Complex spectra
	C, Co                      [][]complex128 // Channels
}

func warperNew(nbuf, osamp, olap, nch int, nanowarp *Nanowarp) (n *warper) {
	// FIXME Only 2x oversampling works, no more, no less.
	nfft := nextpow2(nbuf * osamp)
	nbins := nfft/2 + 1
	n = &warper{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		olap:  olap,
		root:  nanowarp,
	}
	a := &n.a

	makeslices(a, nbins, nfft, nch)

	// Exceptions
	a.Phase = make([]float64, nbins)
	n.arm = make([]bool, nbins)

	blackmanHarris(a.W[:nbuf])

	windowDx(a.W[:nbuf], a.Wd[:nbuf])
	windowT(a.W[:nbuf], a.Wt[:nbuf])
	copy(a.Wr[:nbuf], a.W[:nbuf])
	slices.Reverse(a.Wr[:nbuf])
	n.wgain = windowGain(n.a.W)
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * n.wgain

	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *warper) process3(lin, rin, lout, rout []float64, coeffs, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)

	input := make2(2, len(lin))
	grainbuf := make2(2, n.nfft)
	ingrain := make2(2, n.nfft)
	copy(input[0], lin)
	copy(input[1], rin)
	lastone := 0
	for j := -n.nbuf; ; j += n.hop {
		if j > len(lout)-n.nbuf {
			break
		}
		i := int(phasor[max(0, j)] - float64(n.nbuf/2))

		if i > len(lin)-n.nbuf {
			break
		}
		for ch := range grainbuf {
			clear(ingrain[ch])
			copy(ingrain[ch][max(0, -i):], input[ch][max(0, i):clamp(0, len(lin), i+n.nbuf)])
		}

		c := 1.
		if j > 0 {
			c = 1 / coeffs[j]
			if c != c || math.IsInf(c, 0) {
				c = 1
			}
		}

		n.advance2(ingrain, grainbuf, abs(c), n.root.opts.Quality >= 0 && c == 1)

		// Cut pre-echo in transient regions.
		// TODO Probably won't need this after non-causal PGHI is implemented.
		d := j - lastone
		if c == 1 {
			lastone = j
		} else if d < n.nbuf/2 {
			for ch := range grainbuf {
				rr := grainbuf[ch][max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grainbuf[ch][:n.nbuf/2-d-n.hop], 0)
			}
		}

		if n.root.opts.Onsets && c != 1 {
			for ch := range grainbuf {
				clear(grainbuf[ch])
			}
		}

		loutgrain := lout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		routgrain := rout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]

		add(loutgrain, grainbuf[0][clamp(0, n.nbuf, -j):])
		add(routgrain, grainbuf[1][clamp(0, n.nbuf, -j):])
	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(present, output [][]float64, stretch float64, reset bool) {
	a := &n.a
	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}
	defft := func(out []float64, x []complex128, w bool) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		if w {
			mul(a.S, a.Wr)
		}
		copy(out, a.S)
	}

	clear(a.Mid)
	for ch := range present {
		add(a.Mid, present[ch])
		enfft(a.C[ch], a.W, present[ch])
	}

	enfft(a.X, a.W, a.Mid)
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	// TODO Works ONLY with 2x FFT oversampling.
	olap := float64(n.nbuf / n.hop)
	osampc := float64(n.nfft) / float64(n.nbuf)
	tadv := gettadv(a.X[:n.nbins], a.Xd, olap*osampc)
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch*2./osampc)

	// Forced phase reset.
	bash := func() {
		for w := range a.Phase {
			a.P[w], a.Pphase[w] = cmplx.Polar(a.X[w])
		}
		for ch := range present {
			// This is the only known way to correctly scale gains.
			enfft(a.Co[ch], a.W, present[ch])
			defft(output[ch], a.Co[ch], true)
			// And this works with -40 dB difference.
			// copy(output[ch], present[ch])
			// mul(output[ch], a.W)
			// mul(output[ch], a.Wr)
			// scale(output[ch], (1+n.wgain*float64(n.nbuf)/float64(n.hop))/2)
			// Where this only must be enough:
			// scale(output[ch], n.wgain*float64(n.nbuf)/float64(n.hop))
		}
	}
	if reset {
		bash()
		return
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
		a.E[w] = p
	}
	for ch := range present {
		for w := range a.X {
			a.C[ch][w] /= a.E[w]
		}
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
		a.E[w] = cmplx.Rect(1, a.Phase[w])
	}
	for w := range a.Phase {
		a.Pphase[w] = princarg(a.Phase[w])
	}
	for ch := range output {
		for w := range a.X {
			a.Co[ch][w] = a.C[ch][w] * a.E[w]
		}
		defft(output[ch], a.Co[ch], true)
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
