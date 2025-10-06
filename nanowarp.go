package nanowarp

import (
	"math"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/dsp/fourier"
)

type bufs struct {
	Fr       []float64 // Frame buffer
	S        []float64 // Scratch buffer
	Cs0, Cs1 []complex128

	W, Wd, Wt, Wdt []float64    // Window functions
	X, Xd, Xt, Xdt []complex128 // Complex spectra

	Speckle []float64
	P, Xh   []complex128
}

type Nanowarp struct {
	nfft   int
	hop    int
	outhop int

	q  []float64
	qi int

	fft *fourier.FFT

	norm float64

	a bufs
}

func New() (n *Nanowarp) {
	nfft := 4096
	n = &Nanowarp{
		nfft: nfft,
		hop:  nfft / 32,
	}

	n.a.Fr = make([]float64, nfft)
	n.a.S = make([]float64, nfft)
	n.a.Cs0 = make([]complex128, nfft/2+1)
	n.a.Cs1 = make([]complex128, nfft/2+1)
	n.a.W = make([]float64, nfft)
	n.a.Wd = make([]float64, nfft)
	n.a.Wt = make([]float64, nfft)
	n.a.Wdt = make([]float64, nfft)
	n.a.X = make([]complex128, nfft/2+1)
	n.a.Xd = make([]complex128, nfft/2+1)
	n.a.Xt = make([]complex128, nfft/2+1)
	n.a.Xdt = make([]complex128, nfft/2+1)
	n.a.Speckle = make([]float64, nfft)
	n.a.P = make([]complex128, nfft/2+1)
	n.a.Xh = make([]complex128, nfft/2+1)

	// rn := reflect.ValueOf(n.a)
	// for i := 0; i < rn.NumField(); i++ {
	// 	f := rn.Field(i)
	// 	if f.Kind() == reflect.Slice {
	// 		sz := nfft
	// 		f.
	// 		if f.Index(0).CanComplex() {
	// 			sz = nfft/2 + 1
	// 		}
	// 		f.Grow(sz)
	// 		f.SetLen(sz)
	// 	}
	// }
	n.q = make([]float64, 2*nfft) // Queue is an exception.

	blackman(n.a.W)
	blackmanDx(n.a.Wd)
	windowT(n.a.W, n.a.Wt)
	windowT(n.a.Wd, n.a.Wdt)
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * windowGain(n.a.W)

	n.fft = fourier.NewFFT(nfft)

	return
}

// func (n *Nanowarp) Push(in, out []float64) {
// 	for {
// 		n.fr = append(n.fr, in[:min(len(in), n.nfft-len(n.fr))]...)
// 		in = in[min(len(in), n.hop):]
// 		if len(n.fr) == n.nfft {
// 			n.process()
// 			add(n.q[n.qi:], n.fr)
// 			n.fr = n.fr[:0]
// 			n.qi += n.hop
// 		}
// 		if n.qi > n.nfft {
// 			n.qi = 0
// 			co := min(len(out), n.nfft)
// 			if len(out) > 0 {
// 				copy(out, n.q[:co])
// 				out = out[co:]
// 			} else {
// 				return
// 			}
// 			// //\_____ -> \_______
// 			copy(n.q[0:], n.q[co:])
// 			zero(n.q[n.nfft:])
// 		}
// 	}
// }
/*

p = read1 * (prev/mag(prev))
out = read2 * (p/mag(p))
prev = out

*/

func (n *Nanowarp) Process(in []float64, out []float64, stretch float64) {
	adv := 0
	o := int(math.Floor(float64(n.hop) * stretch))
	for i := o; i < len(in); i += n.hop {
		enfft := func(i int, x []complex128, w []float64) {
			copy(n.a.S, in[i:i+n.nfft])
			mul(n.a.S, w)
			n.fft.Coefficients(x, n.a.S)
		}

		enfft(i, n.a.X, n.a.W)
		enfft(i-o, n.a.Xh, n.a.W)
		enfft(i, n.a.Xd, n.a.Wd)
		enfft(i, n.a.Xt, n.a.Wt)
		enfft(i, n.a.Xdt, n.a.Wdt)

		// Extract per-bin phase reset points through reassignment.
		// Values around zero (after normalization by nfft) are
		// correlated to transient components of the signal.
		speckle(n.a.X, n.a.Xd, n.a.Xt, n.a.Xdt, n.a.Speckle)
		for i := range n.a.Speckle {
			tol1 := min(0.3, 0.17*math.Sqrt(stretch))
			tol2 := min(0.1, 0.05*math.Sqrt(stretch))
			base := float64(n.nfft)
			n.a.Speckle[i] = boolnum(
				n.a.Speckle[i] > base*(1-tol1) && n.a.Speckle[i] < base*(1+tol1) ||
					n.a.Speckle[i] > base*-tol2 && n.a.Speckle[i] < base*+tol2)
			// n.a.Speckle[i] = boolnum(n.a.Speckle[i] < 0)
		}

		// Perform phase vocoding.
		for i := range n.a.Cs0 {
			n.a.Cs0[i] = n.a.P[i] / (n.a.Xh[i] + 1e-15)
		}
		phaselock(n.a.Cs0, n.a.Cs1)
		for i := range n.a.P {
			n.a.P[i] = norm(n.a.Cs0[i]+1e-15) * n.a.X[i]
			if n.a.Speckle[i] > 0.5 {
				n.a.P[i] = n.a.Xh[i]
			}
			n.a.P[i] /= 3
		}

		n.fft.Sequence(n.a.S, n.a.P)
		mul(n.a.S, n.a.W)
		for i := range n.a.S {
			n.a.S[i] /= n.norm
			n.a.S[i] *= stretch
		}
		add(out[adv:adv+n.nfft], n.a.S)
		adv += o
	}
}

// speckle calculates the mixed partial derivative of spectral phase by time and
// frequency from reassigned windowed complex spectra.
// See Fitz, Kelly R., and Sean A. Fulop. “A unified theory of time-frequency reassignment”.
// https://arxiv.org/pdf/0903.3080
func speckle(x, xd, xt, xdt []complex128, out []float64) {
	for i := range x {
		out[i] = real(xdt[i]/x[i]) - real((xt[i]*xd[i])/(x[i]*x[i]))
	}
}

// freqcorr calculates the partial derivative of phase by time.
// Also known as “reassigned frequency correction” or “phase advance”.
func freqcorr(x, xd []complex128, out []float64) {
	for i := range x {
		out[i] = imag(xd[i]/x[i]) / math.Pi
	}
}

// phaselock locks the neighboring phase components of the complex spectrum.
// See Puckette, Miller. “Phase-locked vocoder”
// https://msp.ucsd.edu/Publications/mohonk95.pdf
func phaselock(x []complex128, out []complex128) {
	for i := range x {
		out[i] = x[i] - x[max(0, i-1)] - x[min(len(x)-1, i+1)]
	}
}

func blackman(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = .42 - .5*math.Cos(math.Pi*x*2) + .08*math.Cos(math.Pi*x*4)
	}
}

func blackmanDx(out []float64) {
	for i := range out {
		x := float64(i)/float64(len(out)) + .5
		out[i] = -math.Pi*math.Sin(2*math.Pi*x) - .32*math.Pi*math.Sin(4*math.Pi*x)
	}
}

func windowGain(w []float64) (a float64) {
	for _, e := range w {
		a += e * e
	}
	a /= float64(len(w))
	return
}

func windowT(w, out []float64) {
	n := float64(len(w))
	for i := range w {
		out[i] = w[i] * mix(-n/2, n/2+1, float64(i)/n)
	}
}

func zero[T any](dst []T) {
	var z T
	for i := range dst {
		dst[i] = z
	}
}

func add[T constraints.Float](dst, src []T) {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] += src[i]
	}
}

func mul[T constraints.Float](dst, src []T) {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] *= src[i]
	}
}

func boolnum(b bool) float64 {
	return map[bool]float64{false: 0, true: 1}[b]
}

func mix[F constraints.Float](x, a, b F) F {
	return a*(1-x) + b*x
}

func min[T constraints.Ordered](a ...T) T {
	b := a[0]
	for _, e := range a {
		if e < b {
			b = e
		}
	}
	return b
}

func mag(c complex128) float64 {
	return math.Sqrt(real(c)*real(c) + imag(c)*imag(c))
}

func norm(c complex128) complex128 {
	mag := math.Sqrt(real(c)*real(c) + imag(c)*imag(c))
	return complex(real(c)/mag, imag(c)/mag)
}

func max[T constraints.Ordered](a ...T) T {
	b := a[0]
	for _, e := range a {
		if e > b {
			b = e
		}
	}
	return b
}
