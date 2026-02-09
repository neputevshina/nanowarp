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
)

const (
	minTransientMs = 15
	maxTransientMs = 55
)

type Nanowarp struct {
	fs          int
	left, right []float64

	*warper
	// hpss  *splitter
	*detector

	opts      Options
	stretch   float64
	semitones float64
}

type Options struct {
	// Use time-phase derivative differences from mono in stereo preservation
	// instead of phase differences.
	Diffadv bool
	// Output scaled onsets only.
	Onsets bool
	// Scale the time between transients quadratically.
	Riddim bool
	// Don't perform transient separation, output raw PVDR with phase reset
	// with arbitrary periodicity.
	Raw bool
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

	n.warper = warperNew(4096*w, n) // TODO 6144@48k prob the best

	// a := sargs{}
	// a.nfft = 1024
	// a.hn = 81
	// a.vn = 31
	// a.dilaten = 3
	// a.hq = 0.5
	// a.vq = 0.5
	// a.thresh = 0.75

	// n.hpss = splitterNew(a, float64(int(1)<<w))
	n.detector = detectorNew(1024, samplerate)

	return
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")
	fmt.Fprintln(os.Stderr, "Coeff =", stretch)

	lpfile := make([]float64, len(lin))
	// rpfile := make([]float64, len(lin))

	phasor := make([]float64, len(lout))
	if n.opts.Raw {
		for j := range phasor {
			phasor[j] = float64(j) / stretch
		}
	} else {
		// n.hpss.process(lin, rin, lpfile, rpfile, nil, nil, nil, nil)
		n.detector.process(lin, rin, lpfile)

		fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: conversion")
		n.getPhasor(phasor, lpfile, stretch)
		for j := range phasor {
			phasor[j] = phasor[j] - float64(n.detector.nbuf*2)
		}
	}

	n.warper.process3(lin, rin, lout, rout, phasor, 0)
}

func (n *Nanowarp) getPhasor(phasor, lpfile []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).getPhasor")

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
			for ; i < len(lpfile) && j < len(phasor) && abs(lpfile[i]) > 0; i++ {
				phasor[j] = phasor[j-1] + 1
				j++
			}
			lock = false
		}

	}

	for j := len(phasor) - 1; j >= 0; j-- {
		if phasor[j] >= float64(maxTransientMs*n.fs/1000) {
			phasor[j] = 0
		}
		if phasor[j] >= float64(minTransientMs*n.fs/1000) {
			j -= int(phasor[j])
		}
		if j < 0 {
			break
		}
		phasor[j] = 0
	}

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
				if n.opts.Riddim {
					u = math.Sqrt(u)
				}
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

// func (n *Nanowarp) Push64(l []float64, r []float64) (nl, nr int) {
// 	c := n.nmustcollect()
// 	diffa := func(a, b *[]float64) {
// 		d := len(*a) - c
// 		if d < 0 {
// 			*a = append(*a, (*b)[:min(len(*b), -d)]...)
// 		}
// 	}
// 	diffa(&n.left, &l)
// 	diffa(&n.right, &r)
// 	return 0, 0
// }

// func (n *Nanowarp) Ready() bool {
// 	if len(n.left) != len(n.right) {
// 		panic(`nanowarp: unreachable, stereo buffer length mismatch`)
// 	}
// 	return len(n.left) == n.nmustcollect()
// }

// func (n *Nanowarp) Pull64(left []float64, right []float64) (nl, nr int) {
// 	return 0, 0
// }
