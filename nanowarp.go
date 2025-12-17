package nanowarp

// TODO
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
// + Niemitalo asymmetric windowing?
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
	"fmt"
	"math"
	"os"
	"sync"
)

type Nanowarp struct {
	left, right []float64
	mid         *nanowarp
	side        *nanowarp

	stretch   float64
	semitones float64
}

type nanowarp struct {
	file         []float64
	hfile, pfile []float64
	lower, upper *warper
	hpss         *splitter

	root *Nanowarp
}

const chu = false
const cha = true

type Options struct {
	Masking bool
	Smooth  bool
}

func New(samplerate int, opts Options) (n *Nanowarp) {
	n = &Nanowarp{mid: new(samplerate, &opts), side: new(samplerate, &opts)}
	n.mid.root = n
	n.side.root = n

	return
}

func (n *Nanowarp) Process(in []float64, out []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")
	n.mid.Process(in, out, stretch)
}

func new(samplerate int, opts *Options) (n *nanowarp) {
	n = &nanowarp{}
	// TODO Fixed absolute bandwidth through zero-padding.
	// Hint: nbuf is already there.
	// TODO Find optimal bandwidths.
	w := int(math.Ceil(float64(samplerate) / 48000))
	n.lower = warperNew(4096 * w) // 8192 (4096) @ 48000 Hz // TODO 6144@48k prob the best
	n.lower.masking = opts.Masking
	n.upper = warperNew(64 * w) // 128 (64) @ 48000 Hz
	n.upper.masking = opts.Masking
	// n.hpss = splitterNew(1<<(9+w), float64(int(1)<<w)) // TODO Find optimal size
	n.hpss = splitterNew(512, float64(int(1)<<w), opts.Smooth) // TODO Find optimal size
	return
}

func (n *nanowarp) highestdelay() int {
	return max(n.lower.nbuf, n.upper.nbuf, n.hpss.nbuf)
}

func (n *nanowarp) Process(in []float64, out []float64, stretch float64) {
	n.hfile = make([]float64, len(in))
	n.pfile = make([]float64, len(in))
	n.hpss.process(in, n.pfile, n.hfile)

	wg := sync.WaitGroup{}
	wg.Add(2)
	go func() {
		// TODO Is this delay value correct?
		dc := float64(n.lower.hop-n.upper.hop) * (stretch - 1) * 2
		n.lower.process(n.hfile, out, stretch, dc)
		wg.Done()
	}()
	go func() {
		n.upper.process(n.pfile, out, stretch, 0)
		wg.Done()
	}()
	wg.Wait()
}

/*
nw.Push(left []float32, right []float32) (nl, nr int)
nw.Push64(left []float64, right []float64) (nl, nr int)
nw.PushStereo(buf [][2]float32) (n int)
nw.PushStereo64(buf [][2]float64) (n int)

nw.Pull(left []float32, right []float32) (nl, nr int)
nw.Pull64(left []float64, right []float64) (nl, nr int)
nw.PullStereo(buf [][2]float32) (n int)
nw.PullStereo64(buf [][2]float64) (n int)

nw.Ready() bool
nw.SetTime(stretch float64)
nw.SetPitch(semitones float64)
New() *Nanowarp
*/

func (n *Nanowarp) Push64(left []float64, right []float64) (nl, nr int) {

	return 0, 0
}
