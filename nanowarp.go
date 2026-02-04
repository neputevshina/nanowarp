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
	fs          int
	left, right []float64

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
	n.fs = samplerate
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

const (
	minimumSliceMs          = 50
	minimumSliceStretchedMs = 100
	sliceLahMs              = 30
)

func onsetarg(fs int, onsets []float64) (sa []int) {
	cold := 0
	// prev := 0
	for i, x := range onsets {
		if x > 0 {
			interval := fs * minimumSliceMs / 1000
			if cold < 0 {
				sa = append(sa, i)
				// prev = i
				cold = interval
			}
		}
		cold--
	}
	return
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")

	lhfile := make([]float64, len(lin))
	lpfile := make([]float64, len(lin))
	rhfile := make([]float64, len(lin))
	rpfile := make([]float64, len(lin))
	onsetfile := make([]float64, len(lin))

	n.hpss.process(lin, rin, lpfile, rpfile, lhfile, rhfile, onsetfile, nil)

	if n.opts.Onsets {
		copy(lout, onsetfile)
		copy(rout, onsetfile)
		return
	}

	phasor := make([]int, len(lout))
	fuse := true
	o0 := 0
	o1 := 0
	_ = o1
	_ = o0
	pi := 0
	_ = pi
	for j := range phasor {
		i := int(float64(j) / stretch)
		if lpfile[i] <= 0 && fuse {
			phasor[j] = i
			continue
		}
		fuse = false
		if abs(lpfile[i]) > 0 {
			phasor[j] = phasor[j-1] + 1
		}
	}
	e := 0
	for j := 1; j < len(phasor); j++ {
		if phasor[j] == 0 {
			for k := j; k < len(phasor); k++ {
				if phasor[k] > 0 {
					e = k
					break
				}
			}
			for k := j; k <= e; k++ {
				u := unmix(float64(j), float64(e), float64(k))
				m := mix(float64(phasor[j-1]), float64(e)/stretch, u)
				phasor[k] = int(m)
			}
			j = e
		} else {
			phasor[j] += phasor[e]
		}
	}

	// for i := range lout {
	// 	lout[i] = float64(phasor[i]) / float64(len(lin))
	// 	rout[i] = float64(phasor[i]) / float64(len(lin))
	// }
	// return
	// arg := onsetarg(n.fs, onsetfile)
	// // Shift onsets according to the stretch amount.
	// for i := range arg {
	// 	// arg[i] = max(1, int(float64(arg[i]-(sliceLahMs*n.fs/1000))*stretch))
	// 	arg[i] = int(float64(arg[i]) * stretch)
	// }
	// // Convert onsets to global phasor.
	// phasor := make([]int, len(lout))
	// ai := 1
	// for j := range phasor {
	// 	if j > arg[len(arg)-1] {
	// 		i := int(float64(j) / stretch)
	// 		phasor[j] = i
	// 		continue
	// 	}
	// 	onsa := minimumSliceMs * n.fs / 1000
	// 	tbrk := arg[ai-1] - 1 + onsa
	// 	if j >= arg[ai-1] && j < tbrk {
	// 		phasor[j] = phasor[arg[ai-1]-1] + j - arg[ai-1]
	// 	}
	// 	if j >= tbrk && j <= arg[ai] {
	// 		u := unmix(float64(tbrk), float64(arg[ai]), float64(j))
	// 		m := mix(float64(phasor[tbrk-1]), float64(arg[ai])/stretch, u)
	// 		phasor[j] = int(m)
	// 	}
	// 	if j >= arg[ai] {
	// 		ai++
	// 	}
	// 	// lout[j] = float64(phasor[j]) / float64(len(lin))
	// 	// rout[j] = float64(phasor[j]) / float64(len(lin))
	// }
	// return

	n.lower.process3(lin, rin, lout, rout, phasor, 0)
	return

	wg := sync.WaitGroup{}
	wg.Add(2)
	go func() {
		// dc := float64(n.lower.hop-n.upper.hop) * (stretch - 1) * 2
		// n.lower.process3(lhfile, rhfile, lout, rout, phasor, dc)
		wg.Done()
	}()
	go func() {
		n.upper.process3(lpfile, rpfile, lout, rout, phasor, 0)
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
