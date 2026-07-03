package nanowarp

import (
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"os"
	"slices"

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

func DetectorNew(nfft, fs int, maxTransient, onsetevery int) (n *detector) {
	corr := math.Ceil(float64(fs) / 48000)
	nfft = nfft * int(corr)
	nbuf := nfft
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
	makeslices(&n.a, nbins, nfft, 0, 0)

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

		fill(t, c[1])
		mul(t, n.a.Wr)
		add(ons[i:min(len(lin), i+n.nbuf)], t)
	}

	step := even(int(float64(n.m.maxN) / stretch))
	n.m.Reset(step)
	for i := range ons {
		// Center-windowed dilation
		ons1[max(0, i-step/2)], _ = n.m.Filt(ons[i], bang{})
	}

	for i := range ons1 {
		if ons[i] == ons1[i] {
			onsons = append(onsons, [2]float64{float64(i), ons[i]})
		}
	}

	onsons = append(onsons, [2]float64{float64(len(ons)), 0})

	return
}

func (n *detector) NoveltyCurveProcess(ar dspio.SignalReader, aw dspio.SignalWriter) (err error) {
	fmt.Fprintln(os.Stderr, `(*detector).OnsetFunctionWriter`)

	if gr, ok := ar.(*dspio.GrainReader); ok && gr.Hop != gr.N() {
		panic(`onsetFunctionWriter: non-overlapping reader required`)
	}
	gr := dspio.NewOfflineGrainReader(n.nfft, n.hop, ar)
	gw := dspio.NewOfflineGrainWriter(n.nfft, n.hop, aw)
	gs := make([][]float64, 2)
	for ch := range gs {
		gs[ch] = make([]float64, n.nfft)
	}

	fr := make([]float64, n.nbuf)
	fl := make([]float64, n.nbuf)
	for {
		_, err := gr.SignalRead(nil, gs)
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return err
		}

		c := n.advance(gs[0], gs[1])

		fill(fr, c[0])
		fill(fr[n.hop:], 0)

		fill(fl, c[1])

		// fill(fr, sum(c[:]))

		// mul(fr, n.a.Wr)
		mul(fl, n.a.Wr)
		_, err = gw.SignalWrite(nil, [][]float64{fl, fr})
		if err != nil {
			return err
		}
	}
	// This function is expected to exit when io.EOF is encountered.
}

type onset struct {
	i    int
	pow  float64
	bins int
}

func (n *detector) DilatePeakSelectProcess(ar dspio.SignalReader, aw dspio.SignalWriter, stretch float64, ons chan onset) (err error) {
	fmt.Fprintln(os.Stderr, `(*detector).Dilate`)

	if gr, ok := ar.(*dspio.GrainReader); ok && gr.Hop != gr.N() {
		panic(`onsetFunctionWriter: non-overlapping reader required`)
	}
	step := even(int(float64(n.m.maxN) / stretch))
	gr := dspio.NewOfflineGrainReader(step, step/2, ar)
	gw := dspio.NewOfflineGrainWriter(step, step/2, aw)
	gs := make([][]float64, 2)
	for ch := range gs {
		gs[ch] = make([]float64, step)
	}

	n.m.Reset(step)
	for {
		_, err := gr.SignalRead(nil, gs)
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return err
		}

		for i := range gs[0][:step/2] {
			// Center-windowed dilation
			gs[1][i], _ = n.m.Filt(gs[0][i+step/2], bang{})
		}

		_, err = gw.SignalWrite(nil, [][]float64{gs[1], gs[0]})
		if err != nil {
			return err
		}
	}
	// This function is expected to exit when io.EOF is encountered.
}

func (n *detector) advance(lingrain, ringrain []float64) (activations [4]float64) {
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
			// m := mag(x - px*norm(px/(ppx+1e-10)))
			m := cmplx.Abs(x - px*norm(px*cmplx.Conj(ppx)))
			return m * boolfloat(cmplx.Abs(x) > cmplx.Abs(px))
		}
		a.N[w] = bitsafe(max(
			cnov(a.L[w], a.PL[w], a.PPL[w]),
			cnov(a.R[w], a.PR[w], a.PPR[w])))
	}

	copy(a.PL, a.L)
	copy(a.PR, a.R)
	copy(a.PPL, a.PL)
	copy(a.PPR, a.PR)

	// Crossover frequencies
	// l := hztobin(250, n.nfft, n.fs)
	// m := hztobin(820, n.nfft, n.fs)
	// h := hztobin(2500, n.nfft, n.fs)

	// activations[0] = sum(a.N[:l])
	// activations[1] = sum(a.N[l:m])
	// activations[2] = sum(a.N[m:h])
	// activations[3] = sum(a.N[h:])

	s := sum(a.N)

	// softmax(activations[:])
	// for i, v := range activations {
	// 	activations[i] = bitsafe(v)
	// 	if v < 0.25 {
	// 		activations[i] = 0
	// 	}
	// }

	// activations[0] = float64(
	// 	boolint(activations[0] > 0) |
	// 		boolint(activations[1] > 0)<<1 |
	// 		boolint(activations[2] > 0)<<2 |
	// 		boolint(activations[3] > 0)<<3)
	activations[1] = s

	return
}
