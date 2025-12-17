package nanowarp

import (
	"fmt"
	"math"
	"os"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
)

type splitter struct {
	nfft  int
	nbuf  int
	nbins int
	hop   int
	norm  float64
	corr  float64

	fft  *fourier.FFT
	img  [][]float64
	vimp *mediator[float64, bang]
	himp []*mediator[float64, bang]

	a sbufs
}
type sbufs struct {
	S, Wf, Wr, Wd, Wdt, Wt, M, H, P, A []float64
	Fadv, Pfadv                        []float64
	Xprevs                             [][]complex128 // Lookahead framebuffer
	X, Xd, Xdt, Xt, Y                  []complex128
}

func splitterNew(nfft int, filtcorr float64, smooth bool) (n *splitter) {
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 16

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
		// nhimp := 40 * int(filtcorr)
		nhimp := 21 * int(filtcorr)
		qhimp := 0.75
		n.himp[i] = MediatorNew[float64, bang](nhimp, nhimp, qhimp)
	}
	// nvimp := 21 * int(filtcorr)
	nvimp := 15 * int(filtcorr)
	qvimp := 0.25
	n.vimp = MediatorNew[float64, bang](nvimp, nvimp, qvimp)
	if smooth {
		hann(n.a.Wf)
		copy(n.a.Wr, n.a.Wf)
	} else {
		niemitalo(n.a.Wf)
		// Asymmetric window requires applying reversed copy of itself on synthesis stage.
		copy(n.a.Wr, n.a.Wf)
		slices.Reverse(n.a.Wr)
	}
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

	// floatMatrixToImage(n.img)
	// phasogram := func(name string) {
	// 	e := n.img
	// 	fmt.Println(name)
	// 	if len(e) == 0 {
	// 		fmt.Println(`<skipped>`)
	// 		return
	// 	}
	// 	file, err := os.Create(name)
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	png.Encode(file, floatMatrixToImage(e))
	// }
	// phasogram(fmt.Sprint(rand.Int(), `a.png`))

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

	fadv := func(j int) float64 {
		if mag(a.X[j]) < 1e-6 {
			return 0
		}
		// TODO This phase correction value is guaranteed to be wrong but mostly correct.
		return -real(a.Xt[j]/a.X[j])/float64(n.nbins)*math.Pi - math.Pi/2
	}

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
