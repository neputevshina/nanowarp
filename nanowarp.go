package nanowarp

// TODO
// - Pitch
// - Streaming
// - Time and pitch envelopes
// - Hop size dithering
// - Noise extraction
// + Pre-echo fix
// 	- Niemitalo asymmetric windowing?
//	- See sources of Rubber Band V3
//	- Need dx/dt and dx/dw of it
// - HPSS and lower-upper (req. streaming)
// - Optimizations
//	- Calculate mag(a.X) once
//	- Parallelize through channels

import (
	"container/heap"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"reflect"

	"gonum.org/v1/gonum/dsp/fourier"
)

type bufs struct {
	S, Phase, Pphase []float64    // Scratch buffers
	W, Wd, Wt        []float64    // Window functions
	X, Xd, Xt, P, O  []complex128 // Complex spectra
}

type Nanowarp struct {
	lower, upper *warper
}

func New(samplerate int) (n *Nanowarp) {
	n = &Nanowarp{}
	// TODO Fixed absolute bandwidth through zero-padding.
	// Hint: nbuf is already there.
	// TODO Find optimal bandwidths.
	w := int(math.Ceil(float64(samplerate) / 48000))
	n.lower = warperNew(1 << (13 + w)) // 8192 @ 48000 Hz // TODO 120-150 ms prob the best
	n.upper = warperNew(1 << (9 + w))  // 512 @ 48000 Hz // TODO 12-15 ms
	return
}

func (n *Nanowarp) Process(in []float64, out []float64, stretch float64) {
	n.lower.process(in, out, stretch)
}

type warper struct {
	nfft  int
	nbuf  int
	nbins int
	hop   int

	fft  *fourier.FFT
	arm  []bool
	norm float64
	heap hp

	a bufs
}

func warperNew(nfft int) (n *warper) {
	nbuf := nfft / 2
	nbins := nfft/2 + 1
	olap := 4
	n = &warper{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
	}
	a := &n.a

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
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * windowGain(n.a.W)
	n.arm = make([]bool, nbins)

	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *warper) process(in []float64, out []float64, stretch float64) {
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
			for w := range a.Phase {
				a.Pphase[w] = cmplx.Phase(a.X[w])
			}
			continue
		}

		// Begin of phase gradient heap integration (PGHI).
		a := &n.a
		n.heap = make(hp, n.nbins)
		clear(n.arm)

		for j := range a.X {
			n.arm[j] = true
			// Allow time-phase propagation only for local maxima.
			// Which is a simplest possible auditory masking model.
			// if j == 0 || j == n.nbins-1 ||
			// 	mag(a.X[j-1]) < mag(a.X[j]) && mag(a.X[j+1]) < mag(a.X[j]) {
			n.heap[j] = heaptriple{mag(a.P[j]), j, -1}
		}
		heap.Init(&n.heap)

		// Here we are using time-frequency reassignment[¹] as a way of obtaining
		// phase derivatives. Probably in future these derivatives will be replaced
		// with differences because currently phase is leaking in time domain and
		// finite differences guarrantee idempotency at stretch=1.
		//
		// [¹]: Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
		olap := float64(n.nbuf / n.hop)
		tadv := func(j int) float64 {
			osampc := float64(n.nfft/n.nbuf) / 2
			return (math.Pi*float64(j) + imag(a.Xd[j]/a.X[j])) / (olap * osampc)
		}
		fadv := func(j int) float64 {
			return -real(a.Xt[j]/a.X[j])/float64(n.nfft/2)*math.Pi*stretch - math.Pi/2
		}

		hor := make([]float64, n.nbins)
		ver := make([]float64, n.nbins)
		for len(n.heap) > 0 {
			h := heap.Pop(&n.heap).(heaptriple)
			w := h.w
			switch h.t {
			case -1:
				if n.arm[w] {
					hor[w] = princarg(tadv(w))
					a.Phase[w] = a.Pphase[w] + tadv(w)
					n.arm[w] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w]), w, 0})
				}
			case 0:
				if w > 1 && n.arm[w-1] {
					ver[w-1] = princarg(fadv(w - 1))
					a.Phase[w-1] = a.Phase[w] - fadv(w-1)
					n.arm[w-1] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w-1]), w - 1, 0})
				}
				if w < n.nbins-1 && n.arm[w+1] {
					ver[w+1] = princarg(fadv(w + 1))
					a.Phase[w+1] = a.Phase[w] + fadv(w+1)
					n.arm[w+1] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w+1]), w + 1, 0})
				}
			}
		}

		G[`phasogram.png`] = append(G[`phasogram.png`].([][]float64), hor)
		G[`origphase.png`] = append(G[`origphase.png`].([][]float64), ver)

		for w := range a.Phase {
			m := mag(a.X[w])
			p := a.Phase[w]
			a.O[w] = cmplx.Rect(m, p)
			a.Pphase[w] = princarg(a.Phase[w])
		}
		copy(a.P, a.O)
		// End of PGHI.

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

type heaptriple struct {
	mag  float64
	w, t int
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
