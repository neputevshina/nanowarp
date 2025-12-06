package nanowarp

// TODO
// - Time-domain correctness
// - Pitch
// - Streaming
// - Time and pitch envelopes
// - Hop size dithering
// - Noise extraction
// ± Bubbling artifacts fix
// - Play with HPSS filter sizes and quantiles. Try pre-emphasis maybe?
// - Reset phase accuum sometimes against numerical errors
// - Speed-up is broken
// - Pre-echo fix
// 	- Niemitalo asymmetric windowing?
//	- See sources of Rubber Band V3
//	- Need dx/dt of it
// + HPSS and lower-upper
// - Optimizations
//	+ Calculate mag(a.X) once
//	- Try replacing cmplx.Phase with complex arithmetic and measure the speedup
//	- Replace container.Heap with rankfilt
//	- Parallelize
//		- Parallelize after streaming
//	- Use/port a vectorized FFT library (e.g. SLEEF)
//	- Use float32 (impossible with gonum)
//	- SIMD?

import (
	"container/heap"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"sync"

	"gonum.org/v1/gonum/dsp/fourier"
)

type Nanowarp struct {
	hfile, pfile []float64
	lower, upper *warper
	hpss         *splitter
}

func New(samplerate int) (n *Nanowarp) {
	n = &Nanowarp{}
	// TODO Fixed absolute bandwidth through zero-padding.
	// Hint: nbuf is already there.
	// TODO Find optimal bandwidths.
	w := int(math.Ceil(float64(samplerate)/48000)) - 1
	n.lower = warperNew(1 << (13 + w)) // 8192 (4096) @ 48000 Hz // TODO 6144@48k prob the best
	n.upper = warperNew(1 << (10 + w)) // 1024 (512) @ 48000 Hz // TODO 12-15 ms
	n.hpss = splitterNew(1 << (9 + w)) // TODO Find optimal size
	return
}

func (n *Nanowarp) Process(in []float64, out []float64, stretch float64) {
	n.hfile = make([]float64, len(in))
	n.pfile = make([]float64, len(in))
	n.hpss.process(in, n.pfile, n.hfile)
	// Delay compensation.
	// TODO Streaming.
	// TODO *2 is there as a way to make “bubbling” less noticeable
	dc := n.lower.hop - n.upper.hop
	copy(n.pfile, n.pfile[dc:])

	wg := sync.WaitGroup{}
	wg.Add(2)
	go func() {
		// n.lower.process(n.hfile, out, stretch)
		wg.Done()
	}()
	go func() {
		n.upper.process(n.pfile, out, stretch)
		wg.Done()
	}()
	wg.Wait()

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

	a wbufs
}
type wbufs struct {
	S, M, P, Phase, Pphase []float64    // Scratch buffers
	W, Wd, Wt              []float64    // Window functions
	X, Xd, Xt, O           []complex128 // Complex spectra
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

	makeslices(a, nbins, nfft)

	// Exceptions.
	a.Phase = make([]float64, nbins)
	n.arm = make([]bool, nbins)

	hann(a.W[:nbuf])
	hannDx(a.Wd[:nbuf])
	windowT(a.W[:nbuf], a.Wt[:nbuf])
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * windowGain(n.a.W)

	n.fft = fourier.NewFFT(nfft)

	return
}

// process performs phase gradient heap integration (PGHI) based time stretching.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) process(in []float64, out []float64, stretch float64) {
	a := &n.a
	ih := int(math.Floor(float64(n.hop) / stretch))
	oh := int(math.Floor(float64(n.hop)))
	fmt.Fprintln(os.Stderr, `(*warper).process`)
	fmt.Fprintln(os.Stderr, `inhop:`, ih, `nbuf:`, n.nbuf, `nsampin:`, len(in), `nsampout:`, len(out))
	fmt.Fprintln(os.Stderr, `outhop:`, oh)

	adv := 0
	for i := 0; i < len(in); i += ih {
		enfft := func(i int, x []complex128, w []float64) {
			zero(a.S)
			copy(a.S, in[i:min(len(in), i+n.nbuf)])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)
		enfft(i, a.Xd, a.Wd)
		enfft(i, a.Xt, a.Wt)

		for w := range a.X {
			a.M[w] = mag(a.X[w])
		}
		if i == 0 {
			for w := range a.Phase {
				a.Pphase[w] = cmplx.Phase(a.X[w])
			}
			copy(a.P, a.M)

			for j := range a.S[:n.nbuf] {
				a.S[j] = in[i:min(len(in), i+n.nbuf)][j] * a.W[j] * a.W[j] / n.norm * float64(n.nfft)
			}
			add(out[adv:min(len(out), adv+n.nbuf)], a.S)
			adv += oh

			continue
		}

		a := &n.a
		n.heap = make(hp, n.nbins)
		clear(n.arm)

		for j := range a.X {
			n.arm[j] = true
			// Allow time-phase propagation only for local maxima.
			// Which is a simplest possible auditory masking model.
			if j == 0 || j == n.nbins-1 ||
				a.M[j-1] < a.M[j] && a.M[j+1] < a.M[j] {
				n.heap[j] = heaptriple{a.P[j], j, -1}
			}
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
			return -real(a.Xt[j]/a.X[j])/float64(n.nbins)*math.Pi*stretch - math.Pi/2
		}

		for len(n.heap) > 0 {
			h := heap.Pop(&n.heap).(heaptriple)
			w := h.w
			switch h.t {
			case -1:
				if n.arm[w] {
					a.Phase[w] = a.Pphase[w] + tadv(w)
					n.arm[w] = false
					heap.Push(&n.heap, heaptriple{a.M[w], w, 0})
				}
			case 0:
				if w > 1 && n.arm[w-1] {
					// a.Phase[w-1] = a.Phase[w] - fadv(w-1)
					a.Phase[w-1] = a.Phase[w] + (cmplx.Phase(a.X[w-1])-cmplx.Phase(a.X[w]))*stretch
					_ = fadv
					n.arm[w-1] = false
					heap.Push(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
				}
				if w < n.nbins-1 && n.arm[w+1] {
					// a.Phase[w+1] = a.Phase[w] - fadv(w+1)
					a.Phase[w+1] = a.Phase[w] + (cmplx.Phase(a.X[w+1])-cmplx.Phase(a.X[w]))*stretch
					n.arm[w+1] = false
					heap.Push(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
				}
			}
		}

		copy(a.P, a.M)
		for w := range a.Phase {
			a.O[w] = cmplx.Rect(a.M[w], a.Phase[w])
			a.Pphase[w] = princarg(a.Phase[w])
		}

		n.fft.Sequence(a.S, a.O)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.W)
		add(out[adv:min(len(out), adv+n.nbuf)], a.S)
		adv += oh
	}
}

type splitter struct {
	nfft  int
	nbuf  int
	nbins int
	hop   int
	norm  float64

	fft  *fourier.FFT
	vert *mediator[float64, bang]
	hor  []*mediator[float64, bang]

	a sbufs
}
type sbufs struct {
	S, W, M, H, P []float64
	X, Y          []complex128
}

func splitterNew(nfft int) (n *splitter) {
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 8

	n = &splitter{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
	}
	makeslices(&n.a, nbins, nfft)
	n.hor = make([]*mediator[float64, bang], nbins)

	// TODO Samplerate-independent filter sizes.
	for i := range n.hor {
		nhor := 27
		qhor := 0.5
		n.hor[i] = MediatorNew[float64, bang](nhor, nhor, qhor)
	}
	nvert := 15
	qvert := 0.25
	n.vert = MediatorNew[float64, bang](nvert, nvert, qvert)

	hann(n.a.W)
	n.norm = float64(nfft) * float64(olap) * windowGain(n.a.W)
	n.fft = fourier.NewFFT(nfft)

	return
}

// process performs harmonic-percussive source separation (HPSS).
// See Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.
// (https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)
func (n *splitter) process(in []float64, outp []float64, outh []float64) {
	fmt.Fprintln(os.Stderr, `(*splitter).process`)
	for i := 0; i < len(in); i += n.hop {
		n.vert.Reset(n.vert.N)
		a := &n.a

		enfft := func(i int, x []complex128, w []float64) {
			zero(a.S)
			copy(a.S, in[i:min(len(in), i+n.nbuf)])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)

		for w := range a.X {
			a.M[w] = mag(a.X[w])
		}
		n.vert.filt(a.M, n.vert.N, a.P, REFLECT, 0, 0)
		for w := range a.X {
			m := n.hor[w]
			m.Insert(a.M[w], bang{})
			a.H[w], _ = m.Take()
		}

		for w := range a.X {
			if a.P[w] < a.H[w] {
				a.X[w] = 0
			}
		}

		n.fft.Sequence(a.S, a.X)
		for w := range a.S {
			a.S[w] /= n.norm
		}
		mul(a.S, a.W)
		add(outp[i:min(len(outp), i+n.nbuf)], a.S)

		// Harmonic = original - percussive
		in1 := in[i:min(len(in), i+n.nbuf)]
		for j := range in1 {
			a.S[j] -= in1[j] * a.W[j] * a.W[j] / n.norm * float64(n.nfft)
		}
		add(outh[i:min(len(outh), i+n.nbuf)], a.S)
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
