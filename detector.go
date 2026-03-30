package nanowarp

import (
	"fmt"
	"math"
	"os"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
)

type detector struct {
	nfft        int
	nbuf        int
	nbins       int
	hop         int
	norm, wgain float64
	corr        float64
	thresh      float64
	fs          int

	fft                  *fourier.FFT
	bend, dil, peaks, ms *mediator[float64, bang]
	abox, bbox, cbox     *boxfilt

	a dbufs
}
type dbufs struct {
	S, Wf, Wr, N           []float64
	L, R, PL, PPL, PR, PPR []complex128
}

func detectorNew(nfft, fs sa, maxTransient ms) (n *detector) {
	corr := math.Ceil(float64(fs) / 48000)
	nbuf := nfft * int(corr)
	nbins := nfft/2 + 1
	olap := 16

	n = &detector{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		corr:  corr,
		fs:    fs,
	}
	makeslices(&n.a, nbins, nfft)

	nm := fs * 800 / 1000
	nd := fs * 800 / 1000
	nd2 := fs * 4 * int(maxTransient) / 1000
	nd3 := fs * int(maxTransient) / 1000
	n.bend = mediatorNew[float64, bang](nm, nm, 0.75)
	n.dil = mediatorNew[float64, bang](nd, nd, 1)
	n.peaks = mediatorNew[float64, bang](nd2, nd2, 1)
	n.ms = mediatorNew[float64, bang](nd3, nd3, 1)

	n.abox = boxfiltNew(fs / 10)
	n.bbox = boxfiltNew(fs / 10)
	n.cbox = boxfiltNew(fs / 10)

	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	niemitalo(n.a.Wf)
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	n.wgain = windowGain(n.a.Wf)
	n.norm = float64(nfft) * float64(olap) * n.wgain
	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *detector) process2(lin, rin, ons []float64, stretch float64, onsetevery ms) (onsons [][2]float64) {
	fmt.Fprintln(os.Stderr, `(*detector).process`)

	onsons = make([][2]float64, 0, 1000)

	t := make([]float64, n.nfft)
	for i := 0; i < len(lin); i += n.hop {
		c := n.advance(lin[i:min(len(lin), i+n.nbuf)], rin[i:min(len(lin), i+n.nbuf)])

		fill(t, c)
		mul(t, n.a.Wr)
		add(ons[i:min(len(lin), i+n.nbuf)], t)
	}

	skip := int(onsetevery * ms(n.fs) / 1000 / stretch)
	// skip := int(onsetevery * ms(n.fs) / 1000)
	for i := 0; i < len(ons)-skip; i += skip {
		a := i + argmax(ons[i:i+skip])
		v := ons[a]
		// fill(ons[i:i+skip], 0)
		// ons[a] = v
		onsons = append(onsons, [2]float64{float64(a), v})
	}

	return
}

func (n *detector) advance(lingrain, ringrain []float64) (cnov float64) {
	a := &n.a

	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	enfft(a.L, a.Wf, lingrain)
	enfft(a.R, a.Wf, ringrain)

	for w := range a.L {
		// Complex-domain novelty measure calculation in cartesian form.
		//
		// See Duxbury, C., Bello, J. P., Davies, M., & Sandler, M. (2003, September).
		// Complex domain onset detection for musical signals. In Proc. Digital Audio
		// Effects Workshop (DAFx) (Vol. 1, pp. 6-9). London: Queen Mary University.
		//
		// https://www.audiolabs-erlangen.de/resources/MIR/FMP/C6/C6S1_NoveltyComplex.html
		cnov := func(x, px, ppx complex128) float64 {
			m := mag(x - px*norm(px/(ppx+1e-10)))
			return m * boolfloat(mag(x) > mag(px))
		}
		a.N[w] = bitsafe(max(
			cnov(a.L[w], a.PL[w], a.PPL[w]),
			cnov(a.R[w], a.PR[w], a.PPR[w])))
	}

	copy(a.PL, a.L)
	copy(a.PR, a.R)
	copy(a.PPL, a.PL)
	copy(a.PPR, a.PR)

	cnov = sum(a.N)
	return
}
