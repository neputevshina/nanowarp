package nanowarp

import (
	"container/heap"
	"fmt"
	"math"
	"math/cmplx"
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
	arm  []bool
	norm float64
	heap hp

	a bufs
}

func New() (n *Nanowarp) {
	nfft := 2048
	nbuf := nfft / 2

	nbins := nfft/2 + 1
	olap := 4
	n = &Nanowarp{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		iset:  make(map[int]struct{}, nbins),
		q:     make([]float64, 2*nfft),
	}
	a := &n.a

	a.F = make([][]complex128, 5)

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
	n.arm = make([]bool, nbins)

	n.fft = fourier.NewFFT(nfft)

	return
}

type hp []heaptriple

func (h hp) Len() int           { return len(h) }
func (h hp) Less(i, j int) bool { return h[i].mag > h[j].mag }
func (h hp) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *hp) Push(x any) {
	*h = append(*h, x.(heaptriple))
}
func (h *hp) Pop() any {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func (n *Nanowarp) Process(in []float64, out []float64, stretch float64) {
	a := &n.a
	ih := int(math.Floor(float64(n.hop) / stretch))
	oh := int(math.Floor(float64(n.hop)))
	fmt.Fprintln(os.Stderr, `inhop:`, ih, `nbuf:`, n.nbuf)
	fmt.Fprintln(os.Stderr, `outhop:`, oh)

	adv := 0
	for i := 0; i < len(in); i += ih {
		enfft := func(i int, x []complex128, w []float64) {
			zero(a.S)
			copy(a.S, in[i:i+n.nbuf])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)
		enfft(i, a.Xd, a.Wd)
		enfft(i, a.Xt, a.Wt)

		if i == 0 {
			copy(a.P, a.X)
			continue
		}

		// Begin of PGHI

		a := &n.a
		n.heap = make(hp, n.nbins)
		clear(n.arm)

		for j := range a.X {
			n.arm[j] = true
			n.heap[j] = heaptriple{mag(a.P[j]), j, -1}
		}
		heap.Init(&n.heap)

		olap := float64(n.nbuf / n.hop)
		padv := func(j int) float64 {
			// This phase advance formula is a product of 10 hours of trial and error.
			return (math.Pi*float64(j) + imag(a.Xd[j]/a.X[j])) / olap
		}
		fadv := func(j int) float64 {
			// And this one took two days.
			return -real(a.Xt[j]/a.X[j])/float64(len(a.X))*math.Pi - math.Pi/2
		}

		for len(n.heap) > 0 {
			h := heap.Pop(&n.heap).(heaptriple)
			w := h.w
			switch h.t {
			case -1:
				if n.arm[w] {
					a.Phase[w] = princarg(cmplx.Phase(a.P[w]) + padv(w))
					n.arm[w] = false
					heap.Push(&n.heap, heaptriple{mag(izero(a.X, w)), w, 0})
				}
			case 0:
				if w > 1 && n.arm[w-1] {
					a.Phase[w-1] = princarg(a.Phase[w] - fadv(w))
					n.arm[w-1] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w-1]), w - 1, 0})
				}
				if w < n.nbins-1 && n.arm[w+1] {
					a.Phase[w+1] = princarg(a.Phase[w] + fadv(w))
					n.arm[w+1] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w+1]), w + 1, 0})
				}
			}
		}

		for j := range a.Phase {
			a.A[j] = cmplx.Rect(mag(a.X[j]), a.Phase[j])
		}
		copy(a.P, a.A)

		// End of PGHI

		n.fft.Sequence(a.S, a.P)
		for j := range a.S {
			a.S[j] /= n.norm
			a.S[j] /= stretch
		}
		mul(a.S, a.W)
		add(out[adv:adv+n.nbuf], a.S)
		adv += oh
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

// A stripped-down version of the standard container/heap.

type heaptriple struct {
	mag  float64
	w, t int
}

func less(h *[]heaptriple, i, j int) bool {
	// [...] To build a priority queue, implement the Heap interface with the
	// (negative) priority as the ordering for the Less method, so Push adds
	// items while Pop removes the highest-priority item from the queue.
	return (*h)[i].mag > (*h)[j].mag
}

func heapinit(h *[]heaptriple) {
	// heapify
	n := len(*h)
	for i := n/2 - 1; i >= 0; i-- {
		down(h, i, n)
	}
}

func swap(h *[]heaptriple, i, j int) {
	(*h)[i], (*h)[j] = (*h)[j], (*h)[i]
}

func pop(h *[]heaptriple) heaptriple {
	n := len(*h) - 1
	l := (*h)[n]
	*h = (*h)[:n]
	return l
}

// heappush pushes the element x onto the heap.
// The complexity is O(log n) where n = h.Len().
func heappush(h *[]heaptriple, x heaptriple) {
	*h = append(*h, x)
	up(h, len(*h)-1)
}

// heappop removes and returns the minimum element (according to Less) from the heap.
// The complexity is O(log n) where n = h.Len().
// heappop is equivalent to [heapremove](h, 0).
func heappop(h *[]heaptriple) heaptriple {
	n := len(*h) - 1
	swap(h, 0, n)
	down(h, 0, n)
	return pop(h)
}

// heapremove removes and returns the element at index i from the heap.
// The complexity is O(log n) where n = h.Len().
func heapremove(h *[]heaptriple, i int) any {
	n := len(*h) - 1
	if n != i {
		swap(h, i, n)
		if !down(h, i, n) {
			up(h, i)
		}
	}
	return pop(h)
}

// heapfix re-establishes the heap ordering after the element at index i has changed its value.
// Changing the value of the element at index i and then calling heapfix is equivalent to,
// but less expensive than, calling [heapremove](h, i) followed by a Push of the new value.
// The complexity is O(log n) where n = h.Len().
func heapfix(h *[]heaptriple, i int) {
	if !down(h, i, len(*h)) {
		up(h, i)
	}
}

func up(h *[]heaptriple, j int) {
	for {
		i := (j - 1) / 2 // parent
		if i == j || less(h, j, i) {
			break
		}
		swap(h, i, j)
		j = i
	}
}

func down(h *[]heaptriple, i0, n int) bool {
	i := i0
	for {
		j1 := 2*i + 1
		if j1 >= n || j1 < 0 { // j1 < 0 after int overflow
			break
		}
		j := j1 // left child
		if j2 := j1 + 1; j2 < n && less(h, j2, j1) {
			j = j2 // = 2*i + 2  // right child
		}
		if !less(h, j, i) {
			break
		}
		swap(h, i, j)
		i = j
	}
	return i > i0
}

// G is an implicit global variable map for internal debugging purposes.
// Accesses to this map are either in a thing that is not done yet or safe to remove.
//
// You MUST NOT initialize this map.
// If something fails because of it, it's because Nanowarp is broken,
// and you must use a previous version of the library.
var G map[string]any
