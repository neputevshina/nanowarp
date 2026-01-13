package nanowarp

// TODO
// - Pitch
//	- Good resampler is required
// - Streaming
// - Identical stretch curve for any hop size OR stretch curve sync.
//	- Needed to bring HPSS back for higher quality.
// - Offline dynamic stretch size processing, the last step to high quality
//	- Reset hop near the each transient
//	- Repay time debt after each onset by even parts, do not use exponential decay
//	- Classify material using existing HPSS info by “mostly percussive” and “mostly harmonic”
//		- Simply by summing each sample's amplitude in a windowed grain
//		- Repay more debt on mostly harmonic grains
// - Finite difference scheme for PGHI instead of reassigned
// ~ Time and pitch envelopes
// 	- The internal machinery is already there, just glue pieces together and add UI
// + Hop size dithering
//	+ Harmonic-percussive desync fix (bubbling)
// + Phase drift fix (try long impulse train signal to see it)
// 	+ Time-domain correctness
//	- Probably DTW can help. See https://www.youtube.com/watch?v=JNCVj_RtdZw
//	+ Fixed with transient sync
// - Play with HPSS filter sizes and quantiles. Try pre-emphasis maybe?
// + Reset phase accuum sometimes to counteract numerical errors
//	- Fixed with transient sync
// - Optimizations
//	+ Calculate mag(a.X) once
//	+ Replace container.Heap with rankfilt
//	+ Parallelize
//		- Parallelize in streaming
//	- Use/port a vectorized FFT library (e.g. SLEEF/PFFFT)
//	- Use only float32 (impossible with gonum)
//	- SIMD?
//		- dev.simd branch of Go compiler with intrinsics
//	 - Silent percussive frame elimination

import (
	"fmt"
	"math"
	"os"
	"sync"
)

type Nanowarp struct {
	left, right []float64

	file         []float64
	hfile, pfile []float64
	lower, upper *warper
	hpss         *splitter

	opts      Options
	stretch   float64
	semitones float64
}

type Options struct {
	Diffadv bool
	Onsets  bool
	Finite  bool
}

func New(samplerate int, opts Options) (n *Nanowarp) {
	n = new(samplerate, &opts)
	n.opts = opts

	return
}

func new(samplerate int, opts *Options) (n *Nanowarp) {
	// TODO Fixed absolute bandwidth through zero-padding.
	// Hint: nbuf is already there.
	// TODO Find optimal bandwidths.
	n = &Nanowarp{}
	w := int(math.Ceil(float64(samplerate) / 48000))

	n.lower = warperNew(4096*w, n) // 8192 (4096) @ 48000 Hz // TODO 6144@48k prob the best
	n.lower.diffadv = opts.Diffadv
	n.upper = warperNew(64*w, n) // 128 (64) @ 48000 Hz
	n.lower.diffadv = opts.Diffadv

	n.hpss = splitterNew(512, float64(int(1)<<w))

	return
}

func (n *Nanowarp) nmustcollect() int {
	lop := func(nbuf, hop int) int {
		return nbuf + int(math.Ceil(float64(hop)/n.stretch))
	}
	return max(lop(n.lower.nbuf, n.lower.hop), lop(n.upper.nbuf, n.upper.hop), lop(n.hpss.nbuf, n.hpss.hop))
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")

	lhfile := make([]float64, len(lin))
	lpfile := make([]float64, len(lin))
	rhfile := make([]float64, len(lin))
	rpfile := make([]float64, len(lin))
	onsetfile := make([]float64, len(lin))

	n.hpss.process(lin, rin, lpfile, rpfile, lhfile, rhfile, onsetfile)

	if n.opts.Onsets {
		copy(lout, onsetfile)
		copy(rout, onsetfile)
		return
	}

	wg := sync.WaitGroup{}
	wg.Add(2)
	go func() {
		dc := float64(n.lower.hop-n.upper.hop) * (stretch - 1) * 2
		n.lower.process2(lhfile, rhfile, lout, rout, onsetfile, stretch, dc)
		wg.Done()
	}()
	go func() {
		n.upper.process2(lpfile, rpfile, lout, rout, onsetfile, stretch, 0)
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

func (n *Nanowarp) Push64(l []float64, r []float64) (nl, nr int) {
	c := n.nmustcollect()
	diffa := func(a, b *[]float64) {
		d := len(*a) - c
		if d < 0 {
			*a = append(*a, (*b)[:min(len(*b), -d)]...)
		}
	}
	diffa(&n.left, &l)
	diffa(&n.right, &r)
	return 0, 0
}

func (n *Nanowarp) Ready() bool {
	if len(n.left) != len(n.right) {
		panic(`nanowarp: unreachable, stereo buffer length mismatch`)
	}
	return len(n.left) == n.nmustcollect()
}

func (n *Nanowarp) Pull64(left []float64, right []float64) (nl, nr int) {
	return 0, 0
}
