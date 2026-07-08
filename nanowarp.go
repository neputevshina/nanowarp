package nanowarp

import (
	"fmt"
	"math"
	"os"
	"slices"
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

	// Set algorithm quality.
	//  -2: Don't perform transient separation, output raw PVDR without phase resets.
	//  -1: Extract transients and reset the phase when not stretching.
	//	Introduces clicky artifacts but cleanest for transient-heavy material.
	//	Best numerical stability because of resets.
	//  0:  Same as -1, but detects and bypasses tonal components.
	//	No artifacts, but noticeable slight loss in clarity.
	Quality int

	Hyperparams
}

type Hyperparams struct {
	// Diameter of a transient in milliseconds.
	// Amount of time around the detected transient, for which the signal
	// will be unscaled.
	//
	// Please note that if this value is less than synthesis hop size
	// (1024 samples at 48 kHz sample rate), transient reset is not guaranteed.
	TransientMs int `default:"30"`

	// Size of transient picking filter in milliseconds.
	// Minimum amount of time between two consecutive transient detections.
	PickingMs int `default:"250"`

	// If true, measure pooling size in output time, not in input time.
	//
	// I.e. scale the pooling size with the stretch coefficient.
	// Effective only for stretches, not shrinks, which are always scaled.
	ScalePool bool

	// Minimum amount of frequency bins a trajectory must travel in
	// one frame to be destroyed.
	//
	// Lower values — more vibrato will be detected as transients.
	// On 1 only steady sinusoidal trajectories will be considered tonal.
	HighRidgeHeight int `default:"5"`

	// Minimum amount of time bins a trajectory must travel in total
	// to be considered tonal.
	//
	// Lower values — more tonal preservation and less transient clarity.
	// Higher values — more transient preservation and more interrupts.
	LongRidgeLength int `default:"8"`

	// Maximum radius of influence of each detected tonal trajectory.
	// Phase never be reset at this number of bins around the ridge.
	//
	// Higher values compromise transient quality over tonal quality.
	InfluenceRadius int `default:"3"`

	// The frequency in hertz above which every bin at the transient
	// frame will be reset.
	ResetUpToHz float64 `default:"24000"`
}

func New(samplerate int, opts Options) (n *Nanowarp) {
	structinit(&opts)
	structinit(&opts.Hyperparams)
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
	if opts.PickingMs == 0 {
		opts.PickingMs = 250
	}
	olap := 4

	n.warper = warperNew(4096*w, 2, olap, 2, n)
	n.detector = DetectorNew(1024, samplerate, opts.TransientMs, opts.PickingMs)

	return
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")
	fmt.Fprintln(os.Stderr, "Coeff =", stretch)

	ons := make([]float64, len(lin))
	ons1 := make([]float64, len(lin))

	coeffs := make([]float64, len(lout))
	phasor := make([]float64, len(lout))
	if n.opts.Quality == -2 {
		for j := range coeffs {
			coeffs[j] = 1 / stretch
		}
	} else {
		poolstretch := 1.
		if n.opts.ScalePool || stretch < 1 {
			poolstretch = stretch
		}
		sam := n.detector.process2(lin, rin, ons, ons1, poolstretch)
		// n.getCoeffSignal(coeffs, sam, stretch)

		c := &Curve{}
		n.generatePhasor(c, len(lin), stretch)
		n.bendPhasor(c, sam)

		for j := range phasor {
			fj := float64(j)
			coeffs[j], _ = c.Dy(fj)
			phasor[j], _ = c.ReverseSample(fj)
		}
	}

	n.process6([][]float64{lin, rin}, [][]float64{lout, rout}, coeffs, phasor)
}

func (n *Nanowarp) getCoeffSignal(coeffs []float64, onsets [][2]float64, s float64) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).getCoeffSignal")

	tsa := n.opts.TransientMs * n.fs / 1000

	// Coefficients in returned signal are dt (scan speed, inverse of stretch,
	// which current coefficient describes)

	// Trim transients that are too near to each other
	for k := 0; k < len(onsets)-1; k++ {
		i := onsets[k]
		j := onsets[k+1]
		if j[0]-i[0] < float64(max(n.warper.nbuf, tsa))/s {
			// Leave only louder one
			if i[1] > j[1] {
				onsets[k] = i
			}
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}
	}
	fill(coeffs[:int(onsets[0][0]*s)], 1/s)
	for k := 0; k < len(onsets)-1; k++ {
		i := int(onsets[k][0] * s)
		j := int(onsets[k+1][0] * s)
		if s == 1 {
			fill(coeffs[max(0, i-tsa/2):i+tsa/2], 1)
			fill(coeffs[i+tsa/2:j-tsa/2], 1.00001)
		} else {
			fill(coeffs[max(0, i-tsa/2):i+tsa/2], 1)
			t, x := float64(j-i), float64(tsa)
			r := tsa / 2
			if k == len(onsets) {
				r = 0
			}
			fill(coeffs[i+tsa/2:j-r], (t/s-x)/(t-x))
		}
	}
	fill(coeffs[len(coeffs)-tsa/2:], 1/s)
}

func (n *Nanowarp) generatePhasor(c *Curve, inputLength int, s float64) {
	c.Mutate(func(f []Breakpoint) []Breakpoint {
		in := float64(inputLength)
		return []Breakpoint{Bp(0, 0), Bp(in, in*s)}
	})
}

var println = fmt.Println

func (n *Nanowarp) bendPhasor(c *Curve, onsets [][2]float64) {
	tsa := n.opts.TransientMs * n.fs / 1000
	for k := 0; k < len(onsets)-1; k++ {
		a := onsets[k]
		b := onsets[k+1]
		j, _ := c.Sample(b[0])
		sa, _ := c.Dx(j)
		if b[0]-a[0] < float64(max(n.warper.nbuf, tsa))/sa {
			// Leave only louder one
			if a[1] > b[1] {
				onsets[k] = a
			}
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}
	}
	for k := 0; k < len(onsets)-1; k++ {
		i := onsets[k][0]
		r := float64(tsa / 2)
		j, _ := c.Sample(i)
		a, b := c.Between(i-r), c.Between(i+r)
		c.Mutate(func(f []Breakpoint) []Breakpoint {
			if a != b {
				f = slices.Delete(f, a, b)
			}
			f = slices.Insert(f, a+1, []Breakpoint{{i - r, j - r}, {i + r, j + r}}...)
			return f
		})
	}
}

type Breakpoint struct {
	I, J float64
}

func Bp(i, j float64) Breakpoint { return Breakpoint{I: i, J: j} }

type Curve struct {
	elems       []Breakpoint
	last, rlast int
	start, end  Breakpoint
}

func (c *Curve) Dx(i float64) (v float64, oflow int) {
	if i >= c.end.I {
		return c.dx(len(c.elems) - 2), 1
	}
	if i < c.start.I {
		return c.dx(0), -1
	}
	f := c.Between(i)
	return c.dx(f), 0
}

func (c *Curve) Dy(j float64) (v float64, oflow int) {
	if j >= c.end.J {
		return 1 / c.dx(len(c.elems)-2), 1
	}
	if j < c.start.J {
		return 1 / c.dx(0), -1
	}
	f := c.ReverseBetween(j)
	return 1 / c.dx(f), 0
}

func (c *Curve) dx(f int) float64 {
	delx := (c.elems[f+1].I - c.elems[f].I)
	dely := (c.elems[f+1].J - c.elems[f].J)
	return dely / delx
}

func (c *Curve) Sample(i float64) (j float64, oflow int) {
	if i >= c.end.I {
		return c.end.J, 1
	}
	if i < c.start.I {
		return c.start.J, -1
	}
	f := c.Between(i)
	ni := unmix(c.elems[f].I, c.elems[f+1].I, i)
	j = precisionmix(c.elems[f].J, c.elems[f+1].J, ni)
	return
}

func (c *Curve) ReverseSample(j float64) (i float64, oflow int) {
	if j >= c.end.J {
		return c.end.I, 1
	}
	if j < c.start.J {
		return c.start.I, -1
	}
	f := c.ReverseBetween(j)
	nj := unmix(c.elems[f].J, c.elems[f+1].J, j)
	i = precisionmix(c.elems[f].I, c.elems[f+1].I, nj)
	return
}

func (c *Curve) Between(i float64) (a int) {
	// if c.elems[c.last].I < i {
	// 	c.last = 0
	// }
	if i >= c.end.I {
		return len(c.elems)
	}
	if i < c.start.I {
		return -1
	}
	for f := c.last; f < len(c.elems); f++ {
		if c.elems[f+1].I > i {
			c.last = f
			return f
		}
	}
	panic(`unreachable`)
}

func (c *Curve) ReverseBetween(j float64) (a int) {
	// if c.elems[c.rlast].I < j {
	// 	c.rlast = 0
	// }
	if j >= c.end.J {
		return len(c.elems)
	}
	if j < c.start.J {
		return -1
	}
	for f := c.last; f < len(c.elems); f++ {
		if c.elems[f+1].J > j {
			c.last = f
			return f
		}
	}
	panic(`unreachable`)
}

func (c *Curve) Mutate(f func([]Breakpoint) []Breakpoint) {
	c.elems = f(c.elems)
	c.start = c.elems[0]
	c.end = c.elems[len(c.elems)-1]
	c.last = 0
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
