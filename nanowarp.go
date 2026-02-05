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
	Alt     bool
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

	a := sargs{}
	a.hn = (2*21 - 1)
	a.vn = 15
	if opts.Alt {
		a.hq = 0.5
		a.vq = 0.5
		a.thresh = 0.75
	} else {
		a.hq = 0.75
		a.vq = 0.25
		a.thresh = 0.4
	}
	n.hpss = splitterNew(a, 1024, float64(int(1)<<w))

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

	n.hpss.process(lin, rin, lpfile, rpfile, lhfile, rhfile, onsetfile, nil)

	// ca := 0.
	// j := 0.
	// for i := range lpfile {
	// 	l := abs(lpfile[i])
	// 	r := abs(rpfile[i])
	// 	v := max(l, r)
	// 	if l < 1e-3 || v < ca {
	// 		lpfile[i] = 0
	// 		rpfile[i] = 0
	// 	} else {
	// 		j++
	// 		ca += (v - ca) / float64(j+1)
	// 	}
	// }
	if n.opts.Onsets {
		copy(lout, lpfile)
		copy(rout, rpfile)
		return
	}
	for i := range lpfile {
		lpfile[i] = max(abs(lpfile[i]), abs(rpfile[i]))
	}

	if n.opts.Alt {
		phasor := n.getPhasor(lout, stretch, lpfile)
		for j := range phasor {
			phasor[j] = phasor[j] - float64(n.lower.nbuf/2)
		}
		// for i := 1; i < len(phasor); i++ {
		// 	lout[i] = 1 / (phasor[i] - phasor[i-1])
		// 	if lout[i] < 0 {
		// 		fmt.Println(`!!`, i, lout[i])
		// 		e++
		// 	}
		// 	if lout[i] != lout[i] || math.IsInf(lout[i], 0) {
		// 		fmt.Println(i, lout[i])
		// 		lout[i] = 1
		// 	}
		// }
		// copy(rout, phasor)
		// copy(lout, phasor)
		// println(e)
		n.lower.process3(lin, rin, lout, rout, phasor, 0)
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

const (
	minTransientMs = 15
	maxTransientMs = 50
)

func (n *Nanowarp) getPhasor(lout []float64, stretch float64, lpfile []float64) []float64 {
	phasor := make([]float64, len(lout))
	fuse := true
	o0 := 0.
	lock := true
	for j := 0; j < len(phasor); j++ {
		i := int(float64(j) / stretch)
		if lpfile[i] <= 0 && fuse {
			phasor[j] = float64(j) / stretch
			continue
		}
		fuse = false
		if !lock && abs(lpfile[i]) <= 0 {
			lock = true
		}
		if abs(lpfile[i]) > 0 && lock {
			for ; i < len(lpfile) && abs(lpfile[i]) > 0; i++ {
				phasor[j] = phasor[j-1] + 1
				j++
			}
			lock = false
		}

	}

	// for j := len(phasor) - 1; j >= 0; j-- {
	// 	if phasor[j] >= float64(maxTransientMs*n.fs/1000) {
	// 		phasor[j] = 0
	// 	}
	// 	if phasor[j] >= float64(minTransientMs*n.fs/1000) {
	// 		j -= int(phasor[j])
	// 	}
	// 	phasor[j] = 0
	// }

	for j := range phasor {
		i := int(float64(j) / stretch)
		if phasor[j] > 0 {
			phasor[j] = o0 + phasor[j]
		} else {
			o0 = float64(i)
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
				phasor[k] = m
			}
			if e < j {
				break
			}
			j = e
		}
	}
	for j := e; j < len(phasor); j++ {
		phasor[j] = float64(j) / stretch
	}
	return phasor
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
