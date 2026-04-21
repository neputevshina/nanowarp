package nanowarp

import (
	"fmt"
	"io"
	"math"
	"os"
	"slices"
	"sync/atomic"

	"github.com/neputevshina/nanowarp/dspio"

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
	m           *mediator[float64, bang]

	fft *fourier.FFT

	a dbufs
}
type dbufs struct {
	S, Wf, Wr, N           []float64
	L, R, PL, PPL, PR, PPR []complex128
}

func detectorNew(nfft, fs sa, maxTransient, onsetevery ms) (n *detector) {
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

	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	niemitalo(n.a.Wf)
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	n.wgain = windowGain(n.a.Wf)
	n.norm = float64(nfft) * float64(olap) * n.wgain
	n.fft = fourier.NewFFT(nfft)

	tsa := int(onsetevery) * n.fs / 1000
	n.m = mediatorNew[float64, bang](tsa+1, tsa+1, 1)

	return
}

func (n *detector) process2(lin, rin, ons, ons1 []float64, stretch float64) (onsons [][2]float64) {
	fmt.Fprintln(os.Stderr, `(*detector).process`)

	onsons = make([][2]float64, 0, 1000)

	t := make([]float64, n.nfft)
	for i := 0; i < len(lin); i += n.hop {
		c := n.advance(lin[i:min(len(lin), i+n.nbuf)], rin[i:min(len(lin), i+n.nbuf)])

		fill(t, c)
		mul(t, n.a.Wr)
		add(ons[i:min(len(lin), i+n.nbuf)], t)
	}

	step := even(int(float64(n.m.maxN) / stretch))
	n.m.Reset(step)
	for i := range ons {
		ons1[max(0, i-step/2)], _ = n.m.Filt(ons[i], bang{})
	}

	for i := range ons1 {
		if ons[i] == ons1[i] {
			onsons = append(onsons, [2]float64{float64(i), ons[i]})
		}
	}

	return
}

func (n *detector) onsetFunctionWriter(ar dspio.SignalReader, aw dspio.SignalWriter, stop *atomic.Bool) (err error) {
	fmt.Fprintln(os.Stderr, `(*detector).onsetFunctionWriter`)

	gr := dspio.NewGrainReader(n.nfft, n.hop, ar)
	gw := dspio.NewGrainWriter(n.nfft, n.hop, aw)
	gs := make([][]float64, 2)
	for ch := range gs {
		gs[ch] = make([]float64, n.nfft)
	}

	fr := make([]float64, n.nbuf)
	for {
		if stop.Load() {
			break
		}
		_, err := gr.SignalRead(nil, gs)

		c := n.advance(gs[0], gs[1])

		fill(fr, c)
		mul(fr, n.a.Wr)
		_, err = gw.SignalWrite(nil, [][]float64{fr})
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
	}
	return nil
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
