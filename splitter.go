package nanowarp

import (
	"fmt"
	"math"
	"os"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
)

type splitter struct {
	nfft     int
	nbuf     int
	nbins    int
	hop      int
	norm     float64
	corr     float64
	detector bool

	fft   *fourier.FFT
	img   [][]float64
	vimp  *mediator[float64, bang]
	himp  []*mediator[float64, bang]
	timp1 *mediator[float64, bang]
	timp2 *mediator[float64, bang]
	timp3 *mediator[float64, bang]

	a sbufs
}
type sbufs struct {
	S, Wf, Wr, Wd, Wdt, Wt, M, H, P, A []float64
	Fadv, Pfadv                        []float64
	Xprevs                             [][]complex128 // Lookahead framebuffer
	X, Xd, Xdt, Xt, Y                  []complex128
}

func splitterNew(nfft int, filtcorr float64, _, detector bool) (n *splitter) {
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 16
	fc := int(filtcorr)

	n = &splitter{
		nfft:     nfft,
		nbins:    nbins,
		nbuf:     nbuf,
		hop:      nbuf / olap,
		corr:     filtcorr,
		detector: detector,
	}
	makeslices(&n.a, nbins, nfft)
	n.himp = make([]*mediator[float64, bang], nbins)

	// TODO Log-scale for HPSS and erosion
	for i := range n.himp {
		// nhimp := 40 * fc
		nhimp := 21 * fc
		qhimp := 0.75
		n.himp[i] = MediatorNew[float64, bang](nhimp, nhimp, qhimp)
	}
	// nvimp := 21 * fc
	nvimp := 15 * fc
	qvimp := 0.25
	n.vimp = MediatorNew[float64, bang](nvimp, nvimp, qvimp)

	if detector {
		// TODO Use dedicated and faster dilate filters.
		// Amplitude smoothing filter.
		ntimp1 := 48000 / 50 * fc // Slew wave amplitudes at 50 Hz
		qtimp1 := 0.998           // ≈ Dilate filter
		n.timp1 = MediatorNew[float64, bang](ntimp1, ntimp1, qtimp1)

		// Smoothing quantile filter.
		ntimp2 := 250 * 48 * fc // Release is 250 ms
		qtimp2 := 0.7
		n.timp2 = MediatorNew[float64, bang](ntimp2, ntimp2, qtimp2)

		// Minimum spacing filter.
		// TODO This filter uses only 1s and 0s, optimize appropriately.
		ntimp3 := 20 * 48 * fc
		qtimp3 := 0.99
		n.timp3 = MediatorNew[float64, bang](ntimp3, ntimp3, qtimp3)
	}

	niemitalo(n.a.Wf)
	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	windowT(n.a.Wf, n.a.Wt)
	n.norm = float64(nfft) * float64(olap) * windowGain(n.a.Wf)
	n.fft = fourier.NewFFT(nfft)

	return
}

// process performs harmonic-percussive source separation (HPSS).
// See Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.
// (https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)
func (n *splitter) process(in []float64, pout []float64, hout []float64) {
	fmt.Fprintln(os.Stderr, `(*splitter).process`)
	for i := range n.himp {
		n.himp[i].Reset(n.himp[i].N)
	}

	poutgrain := make([]float64, n.nfft)
	houtgrain := make([]float64, n.nfft)

	for i := 0; i < len(in); i += n.hop {
		n.advanceOld(in[i:min(len(in), i+n.nbuf)], poutgrain, houtgrain)
		add(pout[i:min(len(pout), i+n.nbuf)], poutgrain)
		add(hout[i:min(len(hout), i+n.nbuf)], houtgrain)
	}

	if n.detector {
		prev := 0.
		for i := range in {
			a := math.Abs(pout[i])
			n.timp1.Insert(a, bang{})
			t1, _ := n.timp1.Take()
			n.timp2.Insert(t1, bang{})
			t2, _ := n.timp2.Take()
			a2 := 0.
			if t1 > t2 {
				a2 = 1
			}
			n.timp3.Insert(a2, bang{})
			c, _ := n.timp3.Take()
			if c <= prev {
				pout[i] = 0
			} else {
				pout[i] = 1
			}
			prev = c
		}
	}
}

func (n *splitter) advanceOld(ingrain []float64, poutgrain []float64, houtgrain []float64) {
	n.vimp.Reset(n.vimp.N)
	a := &n.a

	enfft := func(x []complex128, w []float64) {
		clear(a.S)
		copy(a.S, ingrain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	enfft(a.X, a.Wf)

	for w := range a.X {
		a.M[w] = mag(a.X[w])
	}
	n.vimp.filt(a.M, n.vimp.N, a.P, mREFLECT, 0, 0)
	for w := range a.X {
		m := n.himp[w]
		m.Insert(a.M[w], bang{})
		a.H[w], _ = m.Take()
	}

	for w := range a.X {
		if a.P[w] > a.H[w] {
			a.A[w] += 1
		} else {
			a.A[w] = 0
		}
	}

	for w := range a.X {
		if a.A[w] == 0 {
			a.X[w] = 0
		}
	}

	n.fft.Sequence(a.S, a.X)
	for w := range a.S {
		a.S[w] /= n.norm
	}
	mul(a.S, a.Wr)
	copy(poutgrain, a.S)

	// Harmonic = original - percussive
	for j := range ingrain {
		houtgrain[j] = ingrain[j]*a.Wf[j]*a.Wf[j]/n.norm*float64(n.nfft) - a.S[j]
	}
}

func (n *splitter) advanceNew(ingrain []float64, poutgrain []float64, houtgrain []float64) {
	a := &n.a
	enfft := func(x []complex128, w []float64) {
		clear(a.S)
		copy(a.S, ingrain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	enfft(a.X, a.Wf)
	enfft(a.Xt, a.Wt)

	fadv := getfadv(a.X, a.Xt, 1)

	q := 1.0
	for w := 1; w < n.nbins-1; w++ {
		if princarg(fadv(w)) > math.Pi-q &&
			(princarg(fadv(w-1)) > math.Pi-q ||
				princarg(fadv(w+1)) > math.Pi-q) {
			a.A[w-1] = 4
			a.A[w] = 4
			a.A[w+1] = 4
		} else {
			a.A[w] = max(0, a.A[w]-1)
		}
		// if princarg(fadv(w)) > math.Pi-q &&
		// 	princarg(fadv(w-1)) > math.Pi-q &&
		// 	princarg(fadv(w+1)) > math.Pi-q {
		// 	// e[w] = 1
		// } else {
		// 	a.X[w] = 0
		// }
	}
	for w := 0; w < n.nbins; w++ {
		if a.A[w] == 0 || w == 0 {
			a.X[w] = 0
		}
	}

	// n.img = append(n.img, slices.Clone(a.A))

	n.fft.Sequence(a.S, a.X)
	for w := range a.S {
		a.S[w] /= n.norm
	}
	mul(a.S, a.Wr)
	copy(poutgrain, a.S)

	// Harmonic = original - percussive
	for j := range ingrain {
		houtgrain[j] = ingrain[j]*a.Wf[j]*a.Wf[j]/n.norm*float64(n.nfft) - a.S[j]
	}
}
