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
//

import (
	"fmt"
	"math"
	"os"
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
	// Output scaled onsets only.
	Onsets bool

	// Set algorithm quality
	// -1: Don't perform transient separation, output raw PVDR with phase reset
	// each time after arbitrary amount of samples had been elapsed. 4x overlap. Fastest.
	// 0: Extract transients and reset the phase on them. 4x overlap. Slow.
	// 1: Same as 0, but with 8x overlap. Slowest.
	Quality int

	// Time for which signal will be bypassed at the any given transient.
	//
	// If zero will be set to 30.
	TransientMs int

	// The size of the transient picking bucket in milliseconds.
	// An onset is selected at the point of maximum of the novelty function in each bucket.
	//
	// If zero will be set to 300.
	PoolingMs int
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

	if opts.TransientMs == 0 {
		opts.TransientMs = 30
	}
	if opts.PoolingMs == 0 {
		opts.PoolingMs = 300
	}

	// TODO Probably possible to shrink nbuf to 3072 without loss in quality.
	n.warper = warperNew(4096*w, n)
	n.detector = detectorNew(512, samplerate, ms(opts.TransientMs), ms(opts.PoolingMs))

	return
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")
	fmt.Fprintln(os.Stderr, "Coeff =", stretch)

	ons := make([]float64, len(lin))
	ons1 := make([]float64, len(lin))

	coeffs := make([]float64, len(lout))
	phasor := make([]float64, len(lout))
	if n.opts.Quality == -1 {
		for j := range coeffs {
			coeffs[j] = 1 / stretch
		}
	} else {
		sam := n.detector.process2(lin, rin, ons, ons1, stretch)

		// copy(lout, ons)
		// copy(rout, ons1)
		// return

		n.getCoeffSignal(coeffs, sam, stretch)
	}
	for j := range phasor[1:] {
		phasor[j+1] = phasor[j] + coeffs[j+1]
	}

	n.warper.process3(lin, rin, lout, rout, coeffs, phasor)
}

func (n *Nanowarp) getCoeffSignal(coeffs []float64, onsets [][2]float64, s float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).getCoeffs")

	tsa := n.opts.TransientMs * n.fs / 1000

	for k := 0; k < len(onsets)-1; k++ {
		i := onsets[k][0]
		j := onsets[k+1][0]
		if j-i < float64(max(n.warper.nbuf, tsa))/s {
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}
	}
	fill(coeffs[:int(onsets[0][0]*s)], 1/s)
	fill(coeffs[int(onsets[len(onsets)-1][0]*s):], 1/s)
	for k := 0; k < len(onsets)-1; k++ {
		i := int(onsets[k][0] * s)
		j := int(onsets[k+1][0] * s)
		if s == 1 {
			fill(coeffs[max(0, i-tsa/2):i+tsa/2], 1)
			fill(coeffs[i+tsa/2:j-tsa/2], 1.00001)
		} else {
			fill(coeffs[max(0, i-tsa/2):i+tsa/2], 1)
			t, x := float64(j-i), float64(tsa)
			// Coefficients in signal are dt (scan speed, inverse of speedup,
			// which current coefficient describes).
			fill(coeffs[i+tsa/2:j-tsa/2], (t/s-x)/(t-x))
		}

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
