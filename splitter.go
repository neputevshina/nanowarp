package nanowarp

import (
	"fmt"
	"math"
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

	fft   *fourier.FFT
	img   [][]float64
	vimp  *mediator[float64, bang]
	himp  []*mediator[float64, bang]
	timp1 *mediator[float64, bang]
	timp2 *mediator[float64, bang]
	timp3 *mediator[float64, bang]
	tbox  *boxfilt

	a sbufs
}
type sbufs struct {
	S, Wf, Wr, Wd, Wdt, Wt, M, H, P, A []float64
	Fadv, Pfadv                        []float64
	Xprevs                             [][]complex128 // Lookahead framebuffer
	L, R, X, Y                         []complex128
}

func splitterNew(nfft int, filtcorr float64, _, detector bool) (n *splitter) {
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 16
	corr := int(filtcorr)

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
		nhimp := 21 * corr
		qhimp := 0.75
		n.himp[i] = mediatorNew[float64, bang](nhimp, nhimp, qhimp)
	}
	nvimp := 15 * corr
	qvimp := 0.25
	n.vimp = mediatorNew[float64, bang](nvimp, nvimp, qvimp)

	if detector {
		ntimp1 := 50 * 48 * corr
		qtimp1 := 0.5
		n.timp1 = mediatorNew[float64, bang](ntimp1, ntimp1, qtimp1)

		ntimp2 := 50 * 48 * corr
		qtimp2 := 0.98
		n.timp2 = mediatorNew[float64, bang](ntimp2, ntimp2, qtimp2)

		// Minimum spacing filter.
		// TODO This filter uses only 1s and 0s, optimize appropriately.
		ntimp3 := 100 * 48 * corr
		qtimp3 := 0.99
		n.timp3 = mediatorNew[float64, bang](ntimp3, ntimp3, qtimp3)

		n.tbox = boxfiltNew(9 * 48 * corr)
	}

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
}

// extract performs HPSS with transient extraction.
func (n *splitter) extract(lin, rin []float64, lperc, rperc, lharm, rharm, ons []float64) {
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
		n.advanceExtract(lin[i:min(len(lin), i+n.nbuf)], rin[i:min(len(lin), i+n.nbuf)],
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

	n.extractOnsetCurve(ons)
}

func (n *splitter) extractOnsetCurve(ons []float64) {
	ca := 0.
	for i := range ons {
		fi := float64(i)
		if i > 0 {
			ca = (ons[i] + (fi-1)*ca) / fi
		}

		ons[i] -= ca

		n.tbox.Insert(ons[i])
		t1, _ := n.tbox.Take()
		ons[i] = t1 * float64(n.tbox.n)
	}

	for i := range ons {
		if i < 2 {
			continue
		}
		if !(ons[i-1] > ons[i-2] && ons[i-1] > ons[i]) {
			ons[i-2] = 0
		}
		ons[i-2] = max(0, ons[i-2])
	}

	nhop := int(50 * 48 * n.corr)

	for i := 0; i < len(ons); i += nhop {
		win := ons[i:min(len(ons)-1, i+nhop)]
		empty := true
		for j := range win {
			if win[j] > 0 && empty {
				empty = false
				win[j] = 1
			} else {
				win[j] = 0
			}
		}
	}

	// Somehow we get 50 ms of trigger with these parameters.
	const lenphasedropms = 20
	const lahphasedropms = 10

	pdln := int(math.Floor(lenphasedropms * 48 * n.corr))
	pdlah := int(math.Floor(lahphasedropms * 48 * n.corr))

	count := 0
	for i := range ons {
		if ons[clamp(0, len(ons)-1, i+pdlah)] > 0.5 {
			count = pdln
		}
		if count > 0 {
			ons[i] = 0
		} else {
			ons[i] = 1
		}
		count--
	}
}

func (n *splitter) advanceExtract(lingrain, ringrain []float64, lpercgrain, rpercgrain, lharmgrain, rharmgrain, noveltygrain []float64) {
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
			a.A[w] += 1
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
		houtgrain[j] = ingrain[j]*a.Wf[j]*a.Wr[j]/n.norm*float64(n.nfft) - a.S[j]
	}
}
