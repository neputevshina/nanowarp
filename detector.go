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
	nfft  int
	nbuf  int
	nbins int
	hop   int
	fs    int
	m     *mediator[float64, bang]

	fft *fourier.FFT

	a dbufs
}
type dbufs struct {
	S, Wf, Wr              []float64
	N, A, B, X, Y          []float64 `size:"nbins"`
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
		fs:    fs,
	}
	makeslices(&n.a, nbins, nfft, 0, 0)

	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	niemitalo(n.a.Wf)
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	n.fft = fourier.NewFFT(nfft)

	tsa := int(onsetevery) * n.fs / 1000
	n.m = mediatorNew[float64, bang](tsa+1, tsa+1, 1)

	return
}

func (n *detector) process2(lin, rin, ons, ons1 []float64, stretch float64) (onsons []Onset) {
	onsons = make([]Onset, 0, 1000)

	t := make([]float64, n.nfft)
	for i := 0; i < len(lin); i += n.hop {
		c := n.cdodf(lin[i:min(len(lin), i+n.nbuf)], rin[i:min(len(lin), i+n.nbuf)])

		fill(t, c)
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
			onsons = append(onsons, Onset{I: float64(i), Power: ons[i]})
		}
	}

	onsons = append(onsons, Onset{I: float64(len(ons)), Power: 0})

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

		c := n.cdodf(gs[0], gs[1])

		fill(fl, c)
		mul(fl, n.a.Wr)
		_, err = gw.SignalWrite(nil, [][]float64{fl, fr})
		if err != nil {
			return err
		}
	}
	// This function is expected to exit when io.EOF is encountered.
}

type Onset struct {
	I     float64
	Power float64
}

func (n *detector) DilatePeakSelectProcess(ar dspio.SignalReader, aw dspio.SignalWriter, stretch float64, ons chan Onset) (err error) {
	fmt.Fprintln(os.Stderr, `(*detector).Dilate`)

	if gr, ok := ar.(*dspio.GrainReader); ok && gr.Hop != gr.N() {
		panic(`onsetFunctionWriter: non-overlapping reader required`)
	}
	step := even(int(float64(n.m.maxN) / stretch))
	hop := step / 2
	gr := dspio.NewOfflineGrainReader(step, hop, ar)
	var gw *dspio.GrainWriter
	if aw != nil {
		gw = dspio.NewOfflineGrainWriter(step, hop, aw)
	}
	gs := make([][]float64, 2)
	for ch := range gs {
		gs[ch] = make([]float64, step)
	}

	n.m.Reset(step)
	track := 0
	for {
		_, err := gr.SignalRead(nil, gs)
		if err != nil {
			if ons != nil {
				close(ons)
			}
			if err == io.EOF {
				return nil
			}
			return err
		}

		for i := range gs[0][:step/2] {
			// Center-windowed dilation
			gs[1][i], _ = n.m.Filt(gs[0][i+step/2], bang{})
			if gs[1][i] == gs[0][i+step/2] && ons != nil {
				ons <- Onset{I: float64(track + i), Power: gs[1][i]}
			}
		}

		track += hop

		if aw != nil {
			_, err = gw.SignalWrite(nil, [][]float64{gs[1], gs[0]})
			if err != nil {
				return err
			}
		}
	}
	// This function is expected to exit when io.EOF is encountered.
}

// cdodf calculates complex-domain onset detection function for a given stereo grain.
//
// See Duxbury, C., Bello, J. P., Davies, M., & Sandler, M. (2003, September).
// Complex domain onset detection for musical signals. In Proc. Digital Audio
// Effects Workshop (DAFx) (Vol. 1, pp. 6-9). London: Queen Mary University.
//
// https://www.audiolabs-erlangen.de/resources/MIR/FMP/C6/C6S1_NoveltyComplex.html
func (n *detector) cdodf(lingrain, ringrain []float64) (s float64) {
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
		// Cartesian form of CDODF.
		cnov := func(x, px, ppx complex128) float64 {
			m := cmplx.Abs(x - px*norm(px*cmplx.Conj(ppx)))
			return m * boolfloat(cmplx.Abs(x) > cmplx.Abs(px))
		}
		a.N[w] = bitsafe(max(
			cnov(a.L[w], a.PL[w], a.PPL[w]),
			cnov(a.R[w], a.PR[w], a.PPR[w])))
	}

	s = sum(a.N)

	copy(a.PL, a.L)
	copy(a.PR, a.R)
	copy(a.PPL, a.PL)
	copy(a.PPR, a.PR)

	return
}

// superflux calculates an approximation of Superflux onset detection function for a
// given stereo grain.
//
// Not finished.
//
// See Böck, S., & Widmer, G. (2013, September). Maximum filter vibrato suppression
// for onset detection. In Proc. of the 16th Int. Conf. on Digital Audio Effects
// (DAFx). Maynooth, Ireland (Sept 2013) (Vol. 7, p. 4). Citeseer.
//
// https://www.cp.jku.at/research/papers/Boeck_Widmer_DAFx_2013.pdf
func (n *detector) superflux(lingrain, ringrain []float64) (s float64) {
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
		a.X[w] = cmplx.Abs(a.L[w]) + cmplx.Abs(a.R[w])
		a.Y[w] = cmplx.Abs(a.PL[w]) + cmplx.Abs(a.PR[w])
	}

	_ = binfilt(a.X, a.A)
	f := binfilt(a.Y, a.B)

	for n := range f {
		a.N[n] = max(0, a.A[n]-a.B[n])
	}

	s = sum(a.N)

	copy(a.PL, a.L)
	copy(a.PR, a.R)

	return
}

func binfilt(mag, logram []float64) (n int) {
	const scale = 24
	n = int(math.Log2(float64(len(mag)))) * scale
	for i := range n {
		fi := func(i int) float64 { return float64(i) / scale }
		logram[i] = slices.Max(mag[int(math.Pow(2, fi(i))):int(math.Ceil(math.Pow(2, fi(i+1))))])
		logram[i] *= float64(i)
	}
	return
}
