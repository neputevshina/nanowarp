package nanowarp

import (
	"fmt"
	"math"
	"os"
	"reflect"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/dsp/fourier"
)

const eps = 1e-15

type bufs struct {
	Fr                   []float64 // Frame buffer
	S, Phase, S2, S3, S4 []float64 // Scratch buffers
	Cs0, Cs1             []complex128

	W, Wd, Wt, Wdt []float64    // Window functions
	X, Xd, Xt, Xdt []complex128 // Complex spectra
	P, Xn          []complex128
	F              [][]complex128

	Speckle, Pph []float64
	N, A, Pd, Pt []complex128
	Xdh          []complex128
}

type Nanowarp struct {
	nfft        int
	nbuf        int
	nbins       int
	hop         int
	outhop      int
	prevstretch float64
	coeff       float64

	q  []float64
	qi int

	fft  *fourier.FFT
	tol  float64
	iset map[int]struct{}
	norm float64
	heap []heaptriple

	a bufs
}

func New() (n *Nanowarp) {
<<<<<<< HEAD
	nfft := 8192
	nbuf := nfft / 2
=======
	nbuf := 4096
	nfft := nbuf * 2
>>>>>>> master
	nbins := nfft/2 + 1
	n = &Nanowarp{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / 8,
		iset:  make(map[int]struct{}, nbins),
		q:     make([]float64, 2*nfft),
	}
	a := &n.a

	a.F = make([][]complex128, 4)

	// Automatically initialize slices through reflection.
	rn := reflect.ValueOf(&n.a).Elem()
	for i := 0; i < rn.NumField(); i++ {
		f := rn.Field(i)
		if f.Kind() == reflect.Slice {
			c := f.Interface()
			switch c := c.(type) {
			case []complex128:
				f.Set(reflect.ValueOf(make([]complex128, nbins)))
			case []float64:
				f.Set(reflect.ValueOf(make([]float64, nfft)))
			case [][]complex128:
				for i := range c {
					c[i] = make([]complex128, nbins)
				}
			}
		}
	}
	// Exceptions.
	a.Phase = make([]float64, n.nbins)
	hann(a.W[:nbuf])
	hannDx(a.Wd[:nbuf])
	windowT(a.W[:nbuf], a.Wt[:nbuf])
	windowT(a.Wd[:nbuf], a.Wdt[:nbuf])
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * windowGain(n.a.W)

	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *Nanowarp) Process2(in []float64, out []float64, stretch float64) {
	adv := 0
	a := &n.a
	o := int(math.Floor(float64(n.hop) * stretch))
	fmt.Fprintln(os.Stderr, o)
	n.coeff = 5.5

	fmt.Fprintln(os.Stderr, "pghipaper()")
	pghi := n.pghipaper
	for i := 2 * o; i < len(in); i += n.hop {
		enfft := func(i int, x []complex128, w []float64) {
			copy(a.S, in[i:i+n.nfft])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)
		enfft(i, a.Xd, a.Wd)
		enfft(i, a.Xt, a.Wt)

		enfft(i-2*o, a.F[0], a.W)
		enfft(i-1*o, a.F[1], a.W)
		enfft(i, a.F[2], a.W)
		enfft(i+1*o, a.F[3], a.W)

		enfft(i-o, a.P, a.W)
		enfft(i-o, a.Pd, a.Wd)
		enfft(i-o, a.Pt, a.Wt)

		// if i == 2*o {
		// 	for j := range a.Frame {
		// 		a.Frame[j] = angle(a.P[j])
		// 	}
		// }

		// fmt.Fprintln(os.Stderr, "nanowarp.go: phase convergence test")
		// for i := range n.nbins {
		// 	dt := imag(a.Xd[i]/(a.X[i]+eps)) / -math.Pi
		// 	// df := -real(a.Xt[i]/(a.X[i]+eps)) / float64(n.nfft) * 4 * math.Pi
		// 	fmt.Println(dt)
		// }
		// os.Exit(1)
		// enfft(i+o, a.P, a.W)
		// enfft(i+o, a.Pd, a.Wd)
		// enfft(i+o, a.Pt, a.Wt)

		pghi(stretch, a.A)

		n.fft.Sequence(a.S, a.A)
		for i := range a.S {
			a.S[i] /= n.norm
			a.S[i] *= stretch
		}
		mul(a.S, a.W)
		add(out[adv:adv+n.nfft], a.S)
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

func hann(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = 0.5 * (1.0 - math.Cos(2.0*math.Pi*x))
	}
}

func hannDx(out []float64) {
	for i := range out {
		x := float64(i)/float64(len(out)) + .5
		out[i] = math.Pi * math.Sin(2*math.Pi*x)
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

func mix[F constraints.Float](a, b, x F) F {
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

func poltocar(mag, phase float64) complex128 {
	return complex(mag*math.Cos(phase), mag*math.Sin(phase))
}

func csqrt(c complex128) complex128 {
	m := mag(c)
	re := math.Sqrt(0.5 * (m + real(c)))
	im := math.Sqrt(0.5 * (m - real(c)))
	if math.Signbit(imag(c)) {
		im *= -1
	}
	return complex(re, im)
}

func angle(c complex128) float64 {
	return math.Atan2(imag(c), real(c))
}

func norm(c complex128) complex128 {
	mag := math.Sqrt(real(c)*real(c) + imag(c)*imag(c))
	return complex(real(c)/mag, imag(c)/mag)
}

func princarg(phase float64) float64 {
	pi2 := 2 * math.Pi
	return phase - math.Round(phase/pi2)*pi2
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

func izero[T any](sl []T, i int) T {
	var z T
	if i >= len(sl)-1 || i < 0 {
		return z
	}
	return sl[i]
}

func iwrap[T any](sl []T, i int) T {
	if i >= len(sl)-1 || i < 0 {
		i = ((i % len(sl)) + len(sl)) % len(sl)
	}
	return sl[i]
}

// G is an implicit global variable map for internal debugging purposes.
// Accesses to this map are either in a thing that is not done yet or safe to remove.
//
// You MUST NOT initialize this map.
// If something fails because of it, it's because Nanowarp is broken,
// and you must use a previous version of the library.
var G map[string]any
