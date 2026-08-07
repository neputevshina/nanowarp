package nanowarp

import (
	"io"
	"math"
	"math/cmplx"
	"slices"

	"github.com/neputevshina/nanowarp/dspio"

	"gonum.org/v1/gonum/dsp/fourier"
)

type detector struct {
	nfft  int
	nbuf  int
	nch   int
	nbins int
	hop   int
	fs    int
	m     *mediator[float64, bang]

	fft *fourier.FFT

	a dbufs
}
type dbufs struct {
	S, Wf, Wr     []float64
	N, A, B, X, Y []float64 `size:"nbins"`
	L, PL, PPL    [][]complex128
}

func detectorNew(nfft, fs, nch int, onsetevery int) (n *detector) {
	nfft = nextpow2(nfft)
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 16

	n = &detector{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		nch:   nch,
		hop:   nbuf / olap,
		fs:    fs,
	}
	makeslices(&n.a, nbins, nfft, nch, 0)

	// Asymmetric window requires applying reversed copy of itself on synthesis stage.
	niemitalo(n.a.Wf)
	copy(n.a.Wr, n.a.Wf)
	slices.Reverse(n.a.Wr)

	n.fft = fourier.NewFFT(nfft)

	tsa := int(onsetevery) * n.fs / 1000
	n.m = mediatorNew[float64, bang](tsa+1, tsa+1, 1)

	return
}

func (n *detector) noveltyCurveProcess(ar dspio.SignalReader, aw dspio.SignalWriter) (err error) {
	if gr, ok := ar.(*dspio.GrainReader); ok && gr.Hop != gr.N() {
		panic(`onsetFunctionWriter: non-overlapping reader required`)
	}
	gr := dspio.NewOfflineGrainReader(n.nfft, n.hop, ar)
	gw := dspio.NewOfflineOfflineGrainWriter(n.nfft, n.hop, aw)
	gs := make([][]float64, n.nch)
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

		c := n.cdodf(gs)

		fill(fl, c)
		mul(fl, n.a.Wr)
		_, err = gw.SignalWrite(nil, [][]float64{fl, fr})
		if err != nil {
			return err
		}
	}
	// This function is expected to exit when io.EOF is encountered.
}

func (n *detector) dilatePeakSelectProcess(ar dspio.SignalReader, aw dspio.SignalWriter, stretch float64, ons chan Onset) (err error) {
	if gr, ok := ar.(*dspio.GrainReader); ok && gr.Hop != gr.N() {
		panic(`onsetFunctionWriter: non-overlapping reader required`)
	}
	step := even(int(float64(n.m.maxN)))
	hop := step / 2
	gr := dspio.NewOfflineGrainReader(step, hop, ar)
	var gw *dspio.GrainWriter
	if aw != nil {
		gw = dspio.NewOfflineOfflineGrainWriter(step, hop, aw)
	}
	gs := make([][]float64, n.nch)
	for ch := range gs {
		gs[ch] = make([]float64, step)
	}

	n.m.Reset(step)
	track := 0
	cooldown := 0
	if ons != nil {
		defer close(ons)
	}
	for {
		_, err := gr.SignalRead(nil, gs)
		if err != nil {
			// FIXME Last transient in a track can be dropped.
			if err == io.EOF {
				err = nil
			}
			return err
		}

		for i := range gs[0][step/2:] {
			// Center-windowed dilation
			//
			// If process2 (commit 712856d7 and before) was reading ahead like here,
			// the result of that process would be identical to this.
			gs[1][i], _ = n.m.Filt(gs[0][i+step/2], bang{})
			if gs[1][i] == gs[0][i] && track+i-cooldown > step/2 && ons != nil {
				ons <- Onset{I: float64(track + i), Power: gs[1][i]}
				cooldown = track + i
			}
		}

		track += hop

		if aw != nil {
			_, err = gw.SignalWrite(nil, [][]float64{gs[1], gs[0]})
			if err != nil {
				if err == io.EOF {
					err = nil
				}
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
func (n *detector) cdodf(ingrain [][]float64) (s float64) {
	a := &n.a

	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	for ch := range n.nch {
		enfft(a.L[ch], a.Wf, ingrain[ch])
	}

	for w := range a.L[0] {
		// Cartesian form of CDODF.
		cnov := func(x, px, ppx complex128) float64 {
			m := cmplx.Abs(x - px*norm(px*cmplx.Conj(ppx)))
			return m * boolfloat(cmplx.Abs(x) > cmplx.Abs(px))
		}
		nov := 0.
		for ch := range n.nch {
			nov = max(nov, cnov(a.L[ch][w], a.PL[ch][w], a.PPL[ch][w]))
		}
		a.N[w] = bitsafeOrDie(nov)
	}

	s = sum(a.N)

	for ch := range n.nch {
		copy(a.PL[ch], a.L[ch])
		copy(a.PPL[ch], a.PL[ch])
	}

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
// func (n *detector) superflux(lingrain, ringrain []float64) (s float64) {
// 	a := &n.a

// 	enfft := func(x []complex128, w, grain []float64) {
// 		clear(a.S)
// 		copy(a.S, grain)
// 		mul(a.S, w)
// 		n.fft.Coefficients(x, a.S)
// 	}

// 	enfft(a.L, a.Wf, lingrain)
// 	enfft(a.R, a.Wf, ringrain)

// 	for w := range a.L {
// 		a.X[w] = cmplx.Abs(a.L[w]) + cmplx.Abs(a.R[w])
// 		a.Y[w] = cmplx.Abs(a.PL[w]) + cmplx.Abs(a.PR[w])
// 	}

// 	_ = binfilt(a.X, a.A)
// 	f := binfilt(a.Y, a.B)

// 	for n := range f {
// 		a.N[n] = max(0, a.A[n]-a.B[n])
// 	}

// 	s = sum(a.N)

// 	copy(a.PL, a.L)
// 	copy(a.PR, a.R)

// 	return
// }

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
