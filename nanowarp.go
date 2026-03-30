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
	// Time for which signal will be bypassed at the any given transient.
	TransientMs int
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
	n.detector = detectorNew(1024, samplerate, ms(opts.TransientMs))

	return
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")
	fmt.Fprintln(os.Stderr, "Coeff =", stretch)

	ons := make([]float64, len(lin))

	coeffs := make([]float64, len(lout))
	if n.opts.Raw {
		for j := range coeffs {
			coeffs[j] = 1 / stretch
		}
	} else {
		sam := n.detector.process2(lin, rin, ons, stretch, 400)
		n.getCoeffSignal(coeffs, sam, stretch)
	}

	// phs := make([]float64, len(lout))
	// o := coeffs[0]
	// for i := range coeffs[1:] {
	// 	t := coeffs[i+1]
	// 	phs[i+1] = o
	// 	o += t
	// }

	// fill(coeffs, 1.000001)
	// fill(coeffs[3000:4000], 1)

	// copy(lout, ons)
	// copy(lout, phs)
	// println(phs[84391])
	// println(int(phs[len(phs)-1]), len(phs), len(lin))
	// copy(lout, coeffs)
	// return

	n.warper.process3(lin, rin, lout, rout, coeffs, 0)
}

func (n *Nanowarp) getCoeffSignal(coeffs []float64, onsets [][2]float64, s float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).getCoeffs")

	tsa := n.opts.TransientMs * n.fs / 1000

	for k := 0; k < len(onsets)-1; k++ {
		i := onsets[k][0]
		j := onsets[k+1][0]
		if j-i < float64(tsa)/s {
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}
	}
	fill(coeffs[:int(onsets[0][0]*s)], 1/s)
	fill(coeffs[:int(onsets[len(onsets)-1][0]*s)], 1/s)
	for k := 0; k < len(onsets)-1; k++ {
		i := int(onsets[k][0] * s)
		j := int(onsets[k+1][0] * s)
		// i := max(0, int(onsets[k][0]*s)-n.detector.nbuf*2)
		// j := max(0, int(onsets[k+1][0]*s)-n.detector.nbuf*2)
		// i := max(0, int((onsets[k][0]-float64(n.detector.nbuf*2))*s))
		// j := max(0, int((onsets[k+1][0]-float64(n.detector.nbuf*2))*s))
		// i := max(0, int((onsets[k][0]-float64(n.detector.nbuf))*s))
		// j := max(0, int((onsets[k+1][0]-float64(n.detector.nbuf))*s))
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
