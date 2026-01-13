package nanowarp

import (
	"fmt"
	"os"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
)

type splitter struct {
	nfft        int
	nbuf        int
	nbins       int
	hop         int
	norm, wgain float64
	corr        float64
	detector    bool

	fft    *fourier.FFT
	img    [][]float64
	vimp   *mediator[float64, bang]
	himp   []*mediator[float64, bang]
	dilate *mediator[float64, bang]
	tbox   *boxfilt

	a sbufs
}
type sbufs struct {
	S, Wf, Wr, Wd, Wdt, Wt, M, H, P, A []float64
	Fadv, Pfadv                        []float64
	Xprevs                             [][]complex128 // Lookahead framebuffer
	L, R, X, Y                         []complex128
}

func splitterNew(nfft int, filtcorr float64) (n *splitter) {
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 16
	corr := int(filtcorr)

	n = &splitter{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		corr:  filtcorr,
	}
	makeslices(&n.a, nbins, nfft)
	n.himp = make([]*mediator[float64, bang], nbins)

	// TODO Log-scale for HPSS and erosion
	for i := range n.himp {
		nhimp := 21 * corr
		qhimp := 0.75
		n.himp[i] = mediatorNew[float64, bang](nhimp, nhimp, qhimp)
	}
	nvimp := 15 * corr
	qvimp := 0.25
	n.vimp = mediatorNew[float64, bang](nvimp, nvimp, qvimp)

	ntimp2 := 80 * 48 * corr
	qtimp2 := 1.
	n.dilate = mediatorNew[float64, bang](ntimp2, ntimp2, qtimp2)

	n.tbox = boxfiltNew(200 * 48 * corr)

	niemitalo(n.a.Wf)
	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	windowT(n.a.Wf, n.a.Wt)
	n.wgain = windowGain(n.a.Wf)
	n.norm = float64(nfft) * float64(olap) * n.wgain
	n.fft = fourier.NewFFT(nfft)

	return
}

// process performs harmonic-percussive source separation (HPSS) with onset point extraction.
//
// The general method is based on [1], but we are using different quantiles for both horizontal
// and vertical filters. Onsets are extracted by using count of percussive points is a frame as
// a novelty function and then finding peaks using sliding moving maximum in time domain.
//
// [1]: Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.
// (https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)
func (n *splitter) process(lin, rin []float64, lperc, rperc, lharm, rharm, ons []float64) {
	fmt.Fprintln(os.Stderr, `(*splitter).extract`)
	for i := range n.himp {
		n.himp[i].Reset(n.himp[i].N)
	}

	onsgrain := make([]float64, n.nfft)
	nf := n.nfft
	if lperc == nil {
		nf = 0
	}
	lpercgrain := make([]float64, nf)
	rpercgrain := make([]float64, nf)
	lharmgrain := make([]float64, nf)
	rharmgrain := make([]float64, nf)

	for i := 0; i < len(lin); i += n.hop {
		n.advance(lin[i:min(len(lin), i+n.nbuf)], rin[i:min(len(lin), i+n.nbuf)],
			lpercgrain,
			rpercgrain,
			lharmgrain,
			rharmgrain,
			onsgrain)
		if lperc != nil {
			add(lperc[i:min(len(ons), i+n.nbuf)], lpercgrain)
			add(rperc[i:min(len(ons), i+n.nbuf)], rpercgrain)
			add(lharm[i:min(len(ons), i+n.nbuf)], lharmgrain)
			add(rharm[i:min(len(ons), i+n.nbuf)], rharmgrain)
		}
		add(ons[i:min(len(ons), i+n.nbuf)], onsgrain)
	}

	n.onsetCurve(ons)
}

func (n *splitter) onsetCurve(ons []float64) {
	hold := true
	for i := range ons {
		j := i + n.dilate.N/2 // Center the window.
		n.dilate.Insert(ons[min(len(ons)-1, j)], bang{})
		q, _ := n.dilate.Take()
		if ons[i] >= q {
			if hold {
				ons[i] = 1
			} else {
				ons[i] = 0
			}
			hold = false
		} else {
			ons[i] = 0
			hold = true
		}
	}
}

func (n *splitter) advance(lingrain, ringrain []float64, lpercgrain, rpercgrain, lharmgrain, rharmgrain, noveltygrain []float64) {
	n.vimp.Reset(n.vimp.N)
	a := &n.a

	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	enfft(a.L, a.Wf, lingrain)
	enfft(a.R, a.Wf, ringrain)

	for w := range a.X {
		a.M[w] = (mag(a.L[w]) + mag(a.R[w])) / 2
	}
	n.vimp.filt(a.M, n.vimp.N, a.P, mREFLECT, 0, 0)
	for w := range a.X {
		m := n.himp[w]
		m.Insert(a.M[w], bang{})
		a.H[w], _ = m.Take()
	}

	for w := range a.X {
		if a.P[w] > a.H[w] {
			a.A[w] = 1
		} else {
			a.A[w] = 0
		}
	}

	ssum := 0.
	for w := range a.A {
		ssum += a.A[w]
	}
	for w := range a.L {
		if a.A[w] == 0 {
			a.L[w] = 0
			a.R[w] = 0
		}
	}

	for w := range a.S {
		noveltygrain[w] = ssum / n.norm
	}
	mul(noveltygrain, a.Wr)

	if len(lpercgrain) == 0 {
		return
	}

	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.Wr)
		copy(out, a.S)
	}
	defft(lpercgrain, a.L)
	defft(rpercgrain, a.R)

	// Harmonic = original - percussive
	for i := range lingrain {
		w := a.Wf[i] * a.Wr[i] / n.norm * float64(n.nfft)
		lharmgrain[i] = lingrain[i]*w - lpercgrain[i]
		rharmgrain[i] = ringrain[i]*w - rpercgrain[i]
	}

}
