package nanowarp

import (
	"fmt"
	"math"
	"os"
)

// Nanowarp is a high-quality audio time stretching algorithm.
type Nanowarp struct {
	fs          int
	left, right []float64

	*warper
	*detector

	opts      Options
	stretch   float64
	semitones float64
}

// Breakpoint is a phase ramp and coefficient curve breakpoint.
type Breakpoint struct {
	I    int     // Input sample index.
	J    int     // Output sample index.
	X    float64 // Scan speed factor, dJ/dt.
	Type BreakpointType

	fi, fj float64
}

// BreakpointType describes how a breakpoint in the stream
// should be interpolated to the next one.
//
// They correspond to X, J is interpolated as an anti-derivative
// of BreakpointType curve.
type BreakpointType int

const (
	Hold BreakpointType = iota
)

// Fi returns I as float64.
func (b Breakpoint) Fi() float64 { return float64(b.I) }

// Mix advances b to c for the output sample j as the b.BreakpointType curve.
//
// It silently ignores the case when b.J ≰ j ≰ c.J.
func (b Breakpoint) Mix(c Breakpoint, j int) Breakpoint {
	i := 0.
	x := 0.
	switch b.Type {
	case Hold:
		i = unmix(float64(b.J), float64(c.J), float64(j))
		x = 0.
	}
	return Breakpoint{
		I:    int(mix(float64(b.I), float64(c.I), i)),
		J:    j,
		X:    mix(b.X, c.X, x),
		Type: b.Type,
		fi:   mix(b.fi, c.fi, i),
		fj:   mix(b.fj, c.fj, i),
	}
}

type Options struct {
	// Output scaled onsets only.
	Onsets bool

	// Set algorithm quality
	// -1: Don't perform transient separation, output raw PVDR with phase reset
	// every arbitrary amount of seconds. 4x overlap. Fastest.
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
	sam := []Breakpoint{{I: 0, J: 0, X: 1 / stretch, Type: Hold}}
	if n.opts.Quality > -1 {
		sam = n.detector.process2(lin, rin, ons, ons1, stretch)
		n.getCoeffSignal(coeffs, sam, stretch)
		sam = n.convertOnsets(sam, 1/stretch)
	}

	_ = phasor

	// for i := range sam {
	// 	lout[sam[i].J] = sam[i].X
	// }

	// println(sam)

	for j := range phasor[1:] {
		phasor[j+1] = phasor[j] + coeffs[j+1]
	}

	n.warper.process3(lin, rin, lout, rout, coeffs, phasor, sam)
}

func (n *Nanowarp) getCoeffSignal(coeffs []float64, onsets []Breakpoint, s float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).getCoeffs")

	tsa := n.opts.TransientMs * n.fs / 1000

	for k := 0; k < len(onsets)-1; k++ {
		i := onsets[k].I
		j := onsets[k+1].I
		if float64(j-i) < float64(max(n.warper.nbuf, tsa))/s {
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}
	}
	fill(coeffs[:int(onsets[0].Fi()*s)], 1/s)
	fill(coeffs[int(onsets[len(onsets)-1].Fi()*s):], 1/s)
	for k := 0; k < len(onsets)-1; k++ {
		i := int(onsets[k].Fi() * s)
		j := int(onsets[k+1].Fi() * s)
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

func (n *Nanowarp) convertOnsets(onsets []Breakpoint, s float64) (coeff []Breakpoint) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).convertOnsets")

	tsa := n.opts.TransientMs * n.fs / 1000

	add := func(i, j int, x float64, fi, fj float64) {
		if fi < 0 {
			fi = float64(i)
		}
		if fj < 0 {
			fj = float64(j)
		}
		coeff = append(coeff, Breakpoint{I: i, J: j, X: x, fi: fi, fj: fj})
	}
	add(0, 0, s, -1, -1)
	for k := 0; k < len(onsets)-1; k++ {
		if float64(onsets[k+1].I-onsets[k].I) < float64(max(n.warper.nbuf, tsa))*s {
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}

		l := int(onsets[k].Fi() / s)
		r := int(onsets[k+1].Fi() / s)
		t, x := float64(r-l), float64(tsa)
		add(onsets[k].I-tsa/2, l-tsa/2, 1, -1, onsets[k].Fi()/s-x/2)
		add(onsets[k].I+tsa/2, l+tsa/2, (t*s-x)/(t-x), -1, onsets[k].Fi()/s+x/2)
	}
	coeff[0].X = float64(coeff[1].I) / float64(coeff[1].J)
	coeff[len(coeff)-1].X = s
	coeff = append(coeff, coeff[len(coeff)-1])
	return
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
