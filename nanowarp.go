package nanowarp

// TODO
// - Xsynth (nfft=1024) time-scaled with nearest-neighbor stretched magnitude spectrum
//	to preserve tonal balance.
// - Pitch
//	- Good resampler is required
// - Streaming
// - Time and pitch envelopes
// + Replace HPSS with phase-based impulse detection
// - 2048 window, 4096 for bass
// + Hop size dithering
//	+ Harmonic-percussive desync fix (bubbling)
// - Phase drift fix (try long impulse train signal to see it)
// 	+ Time-domain correctness
//	- Probably DTW can help. See https://www.youtube.com/watch?v=JNCVj_RtdZw
// - HPSS erosion
//	- Empty frame elimination
// - Play with HPSS filter sizes and quantiles. Try pre-emphasis maybe?
// - Reset phase accuum sometimes to counteract numerical errors
// - Try to simply resample percussive content by stretch size
//	- See Royer, T. (2019). Pitch-shifting algorithm design and applications in music.
// + WAV float32 input and output
// - Niemitalo asymmetric windowing?
//	- See sources of Rubber Band V3
//	- Need dx/dt of it
// + HPSS and lower-upper
// - Optimizations
//	+ Calculate mag(a.X) once
//	- Replace container.Heap with rankfilt
//	+ Parallelize
//		- Parallelize in streaming
//	- Use/port a vectorized FFT library (e.g. SLEEF)
//	- Use only float32 (impossible with gonum)
//	- SIMD?
//		- dev.simd branch of Go compiler with intrinsics
// - From tests, zplane Elastiqué seems to use 3072-point-in-time window for everything?

import (
	"container/heap"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"math/cmplx"
	"math/rand"
	"os"
	"slices"
	"sync"

	"gonum.org/v1/gonum/dsp/fourier"
)

type Nanowarp struct {
	hfile, pfile []float64
	lower, upper *warper
	hpss         *splitter
}

const chu = true
const cha = false

func New(samplerate int) (n *Nanowarp) {
	n = &Nanowarp{}
	// TODO Fixed absolute bandwidth through zero-padding.
	// Hint: nbuf is already there.
	// TODO Find optimal bandwidths.
	w := int(math.Ceil(float64(samplerate) / 48000))
	n.lower = warperNew(4096 * w) // 8192 (4096) @ 48000 Hz // TODO 6144@48k prob the best
	n.upper = warperNew(64 * w)   // 256 (128) @ 48000 Hz
	// n.hpss = splitterNew(1<<(9+w), float64(int(1)<<w)) // TODO Find optimal size
	n.hpss = splitterNew(2048, float64(int(1)<<w)) // TODO Find optimal size
	return
}

func (n *Nanowarp) Process(in []float64, out []float64, stretch float64) {
	// n.lower.process(in, out, stretch, 0)
	// return
	if stretch >= 1 {
		n.hfile = make([]float64, len(in))
		n.pfile = make([]float64, len(in))
		n.hpss.process(in, n.pfile, n.hfile)

		wg := sync.WaitGroup{}
		wg.Add(2)
		go func() {
			// TODO Is this delay value correct?
			n.lower.process(n.hfile, out, stretch, float64(n.lower.hop-n.upper.hop)*(stretch-1)*2)
			wg.Done()
		}()
		go func() {
			n.upper.process(n.pfile, out, stretch, 0)
			wg.Done()
		}()
		wg.Wait()

	} else {
		// TODO Too lazy, do something more smart
		n.lower.process(in, out, stretch, 0)
	}
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
	img  [][]float64

	a wbufs
}
type wbufs struct {
	S, M, P, Phase, Pphase []float64 // Scratch buffers
	Fadv, Pfadv            []float64
	Xprevs                 []complex128 // Lookahead framebuffer
	W, Wd, Wt              []float64    // Window functions
	X, Xd, Xt, O           []complex128 // Complex spectra
}

func warperNew(nbuf int) (n *warper) {
	nextpow2 := int(math.Floor(math.Pow(2, math.Ceil(math.Log2(float64(nbuf))))))
	nfft := nextpow2 * 2
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

func (n *warper) process(in []float64, out []float64, stretch float64, delay float64) {
	inhop := float64(n.hop) / stretch
	ih, fh := math.Modf(inhop)
	fmt.Fprintln(os.Stderr, `(*warper).process`)
	fmt.Fprintln(os.Stderr, `stretch:`, stretch, `nbuf:`, n.nbuf, `nsampin:`, len(in), `nsampout:`, len(out))
	fmt.Fprintln(os.Stderr, `inhop:`, inhop, `whole:`, ih, `frac:`, fh, `interval:`, float64(n.hop)/(ih+1), `-`, float64(n.hop)/(ih))
	fmt.Fprintln(os.Stderr, `outhop:`, n.hop)

	n.start(in, out)
	outgrain := make([]float64, n.nfft)

	id, fd := math.Modf(delay)
	j := n.hop + int(id)
	dh := fh
	dd := fd
	for i := int(ih); i < len(in); i += int(ih) {
		if j > len(out) {
			break
		}
		dh += fh
		if dh > 0 {
			i += 1
			dh -= 1
			n.advance(in[i:min(len(in), i+n.nbuf)], outgrain, float64(n.hop)/(ih+1))
		} else {
			n.advance(in[i:min(len(in), i+n.nbuf)], outgrain, float64(n.hop)/(ih))
		}
		dd += fd
		if dh > 0 {
			j += 1
			dd -= 1
		}
		add(out[j:min(len(out), j+n.nbuf)], outgrain)
		j += n.hop
	}

	if cha {
		floatMatrixToImage(n.img)
		phasogram := func(name string) {
			e := n.img
			fmt.Println(name)
			if len(e) == 0 {
				fmt.Println(`<skipped>`)
				return
			}
			file, err := os.Create(name)
			if err != nil {
				panic(err)
			}
			png.Encode(file, floatMatrixToImage(e))
		}
		phasogram(fmt.Sprint(rand.Int(), `a.png`))
	}

}

func (n *warper) start(in []float64, out []float64) {
	a := &n.a

	clear(a.S)
	copy(a.S, in[:min(len(in), n.nbuf)])
	mul(a.S, a.W)
	n.fft.Coefficients(a.X, a.S)

	// fadv := getfadv(a.X, a.Xt, 1)
	for w := range a.X {
		a.M[w] = mag(a.X[w])
		a.Pphase[w] = cmplx.Phase(a.X[w])
		// a.Pfadv[w] =
	}
	copy(a.P, a.M)

	for j := range a.S[:n.nbuf] {
		a.S[j] = in[:min(len(in), n.nbuf)][j] * a.W[j] * a.W[j] / n.norm * float64(n.nfft)
	}
	add(out[:min(len(out), n.nbuf)], a.S)
}

func getfadv(x, xt []complex128, stretch float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) < 1e-6 {
			return 0
		}
		// TODO This phase correction value is guaranteed to be wrong but is mostly correct.
		return -real(xt[j]/x[j])/float64(len(x))*math.Pi*stretch - math.Pi/2

	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(ingrain []float64, outgrain []float64, stretch float64) {
	a := &n.a
	enfft := func(x []complex128, w []float64) {
		clear(a.S)
		copy(a.S, ingrain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	enfft(a.X, a.W)
	enfft(a.Xd, a.Wd)
	enfft(a.Xt, a.Wt)

	for w := range a.X {
		a.M[w] = mag(a.X[w])
	}

	n.heap = make(hp, n.nbins)
	clear(n.arm)

	for j := range a.X {
		n.arm[j] = true
		// // Disabled, adds more pre-echo.
		// // Allow time-phase propagation only for local maxima.
		// // Which is a simplest possible auditory masking model.
		// if j == 0 || j == n.nbins-1 ||
		// 	a.M[j-1] < a.M[j] && a.M[j+1] < a.M[j] {
		n.heap[j] = heaptriple{a.P[j], j, -1}
		// }
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
		if mag(a.X[j]) < 1e-6 {
			return 0
		}
		// TODO osampc is wrong. Does not work with non-power of 2 FFT sizes?
		osampc := float64(n.nfft/n.nbuf) / 2
		return (math.Pi*float64(j) + imag(a.Xd[j]/a.X[j])) / (olap * osampc)

	}
	fadv := func(j int) float64 {
		if mag(a.X[j]) < 1e-6 {
			return 0
		}
		// TODO This phase correction value is guaranteed to be wrong but is mostly correct.
		return -real(a.Xt[j]/a.X[j])/float64(n.nbins)*math.Pi*stretch - math.Pi/2

	}

	// e := slices.Repeat([]float64{0}, n.nbins)
	for w := range a.X {
		a.Fadv[w] = princarg(fadv(w))
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
				a.Phase[w-1] = a.Phase[w] - a.Fadv[w-1]
				n.arm[w-1] = false
				heap.Push(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				a.Phase[w+1] = a.Phase[w] + a.Fadv[w+1]
				n.arm[w+1] = false
				heap.Push(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
			}
		}
	}

	// e := make([]float64, n.nbins)
	// Match phase reset points.
	// See “Axel Roebel. A new approach to transient processing in the phase vocoder.
	// 6th International Conference on Digital Audio Effects (DAFx), Sep 2003”
	// for w := 1; w < n.nbins-1; w++ {
	// 	m := func(j int) bool { return a.Pfadv[j]-a.Fadv[j] > math.Pi }
	// 	if m(w) && m(w-1) && m(w+1) {
	// 		a.Phase[w-1] = cmplx.Phase(a.X[w-1])
	// 		a.Phase[w+1] = cmplx.Phase(a.X[w+1])
	// 		a.Phase[w] = cmplx.Phase(a.X[w])
	// 		// e[w] = 1
	// 		// e[w-1] = 1
	// 		// e[w+1] = 1
	// 	}
	// }
	copy(a.Pfadv, a.Fadv)
	// n.img = append(n.img, e)

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

	copy(outgrain, a.S)
}

type splitter struct {
	nfft  int
	nbuf  int
	nbins int
	hop   int
	norm  float64
	corr  float64

	fft *fourier.FFT
	img [][]float64

	a sbufs
}
type sbufs struct {
	S, W, Wd, Wdt, Wt, M, H, P, A []float64
	Fadv, Pfadv                   []float64
	Xprevs                        [][]complex128 // Lookahead framebuffer
	X, Xd, Xdt, Xt, Y             []complex128
}

func splitterNew(nfft int, filtcorr float64) (n *splitter) {
	nbuf := nfft
	nbins := nfft/2 + 1
	olap := 8

	n = &splitter{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		corr:  filtcorr,
	}
	makeslices(&n.a, nbins, nfft)

	hann(n.a.W)
	windowT(n.a.W, n.a.Wt)

	hannDx(n.a.Wd)
	windowT(n.a.Wd, n.a.Wdt)

	n.norm = float64(nfft) * float64(olap) * windowGain(n.a.W)
	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *splitter) process(in []float64, pout []float64, hout []float64) {
	fmt.Fprintln(os.Stderr, `(*splitter).process`)

	poutgrain := make([]float64, n.nfft)
	houtgrain := make([]float64, n.nfft)

	for i := 0; i < len(in); i += n.hop {
		n.advance(in[i:min(len(in), i+n.nbuf)], poutgrain, houtgrain)
		add(pout[i:min(len(pout), i+n.nbuf)], poutgrain)
		add(hout[i:min(len(hout), i+n.nbuf)], houtgrain)
	}
	if chu {
		floatMatrixToImage(n.img)
		phasogram := func(name string) {
			e := n.img
			fmt.Println(name)
			if len(e) == 0 {
				fmt.Println(`<skipped>`)
				return
			}
			file, err := os.Create(name)
			if err != nil {
				panic(err)
			}
			png.Encode(file, floatMatrixToImage(e))
		}
		phasogram(fmt.Sprint(rand.Int(), `a.png`))
	}
}

func (n *splitter) advance(ingrain []float64, poutgrain []float64, houtgrain []float64) {
	a := &n.a
	enfft := func(x []complex128, w []float64) {
		clear(a.S)
		copy(a.S, ingrain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	const nTransientFrames = 4

	enfft(a.X, a.W)
	// Lookahead delay
	if len(a.Xprevs) <= nTransientFrames {
		a.Xprevs = append(a.Xprevs, a.X)
		return
	}
	copy(a.Xprevs, a.Xprevs[1:])
	copy(a.Xprevs[nTransientFrames-1], a.X)

	enfft(a.Xt, a.Wt)
	enfft(a.Xd, a.Wd)
	enfft(a.Xdt, a.Wdt)

	fadv := getfadv(a.X, a.Xt, 1)
	for w := range a.X {
		a.Fadv[w] = princarg(fadv(w))
		// 	if mag(a.X[w]) == 0 {
		// 		a.A[w] = 0
		// 	} else {
		// 		a.A[w] = real(a.Xdt[w]/a.X[w]) - real((a.Xd[w]*a.Xt[w])/(a.X[w]*a.X[w]))
		// 		if a.A[w] > 0 {
		// 			a.A[w] = math.Log(a.A[w])
		// 		} else {
		// 			a.A[w] = math.Log(-a.A[w])
		// 		}
		// 	}
	}

	for w := 1; w < n.nbins-1; w++ {
		// a.A[w] = a.Pfadv[w] - a.Fadv[w]
		m := func(j int) bool { return a.Pfadv[j]-a.Fadv[j] > math.Pi }
		if m(w) && m(w-1) && m(w+1) {
			n := float64(nTransientFrames)
			a.A[w] = n
			a.A[w-1] = n
			a.A[w+1] = n
		} else {
			a.A[w] = max(0, a.A[w]-1)
		}
	}
	copy(a.Pfadv, a.Fadv)

	// for w := 1; w < n.nbins-1; w++ {
	// 	m := func(j int) bool { return a.Pfadv[j]-a.Fadv[j] > math.Pi }
	// 	if m(w) && m(w-1) && m(w+1) {
	// 		n := float64(nTransientFrames)
	// 		a.A[w] = n
	// 		a.A[w-1] = n
	// 		a.A[w+1] = n
	// 	} else {
	// 		a.A[w] = max(0, a.A[w]-1)
	// 	}
	// }
	// copy(a.Pfadv, a.Fadv)

	// q := 1.0
	// for w := 1; w < n.nbins-1; w++ {
	// 	if princarg(a.Fadv[w]) > math.Pi-q &&
	// 		(princarg(a.Fadv[w-1]) > math.Pi-q ||
	// 			princarg(a.Fadv[w+1]) > math.Pi-q) {
	// 		a.A[w-1] = 4
	// 		a.A[w] = 4
	// 		a.A[w+1] = 4
	// 	} else {
	// 		a.A[w] = max(0, a.A[w]-1)
	// 	}
	// }
	for w := 0; w < n.nbins; w++ {
		a.X[w] = a.Xprevs[0][w]
		if a.A[w] == 0 || w == 0 {
			a.X[w] = 0
		}
	}
	if chu {
		n.img = append(n.img, slices.Clone(a.A))
	}

	n.fft.Sequence(a.S, a.X)
	for w := range a.S {
		a.S[w] /= n.norm
	}
	mul(a.S, a.W)
	copy(poutgrain, a.S)

	// Harmonic = original - percussive
	for j := range ingrain {
		houtgrain[j] = ingrain[j]*a.W[j]*a.W[j]/n.norm*float64(n.nfft) - a.S[j]
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

func floatMatrixToImage(data [][]float64) image.Image {
	if len(data) == 0 || len(data[0]) == 0 {
		return nil
	}

	height := len(data[0])
	width := len(data)

	// Find min and max
	minVal := math.Inf(1)
	maxVal := math.Inf(-1)
	for _, row := range data {
		for _, v := range row {
			if v < minVal {
				minVal = v
			}
			if v > maxVal {
				maxVal = v
			}
		}
	}
	fmt.Println(minVal, maxVal)

	scale := 1.
	offset := 3.14
	scale = 255.0 / (maxVal - minVal)
	offset = -minVal * scale

	img := image.NewGray(image.Rect(0, 0, width, height))
	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			val := data[x][y]*scale + offset
			if val < 0 {
				val = 0
			}
			if val > 255 {
				val = 255
			}
			img.SetGray(x, y, color.Gray{Y: uint8(val + 0.5)})
		}
	}

	return img
}
