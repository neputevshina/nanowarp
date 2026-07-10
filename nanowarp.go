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

	warper   *warper
	detector *detector

	opts      Options
	stretch   float64
	semitones float64
}

type Options struct {
	// Output scaled onsets only.
	Onsets bool

	// Set algorithm quality.
	//  -2: Don't perform transient separation, output raw PVDR without phase resets.
	//	Still resets the phase if derivative of input timemap is 1.
	//	Useful if timemap is generated using external onset detector.
	//  -1: Extract transients and reset the phase when not stretching.
	//	Introduces clicky artifacts but cleanest for transient-heavy material.
	//	Best numerical stability because of full-frame resets.
	//  0:  Same as -1, but detects and bypasses tonal components.
	//	No artifacts, but noticeable slight loss in clarity.
	Quality int

	// Channel for receiving processing progress.
	//
	// If not nil, this channel will receive current input and output sample index
	// pair for every 5 seconds of output and at the start and end of processing.
	//
	// Nanowarp will close the channel after the end of processing.
	Progress chan<- Breakpoint

	Hyperparams
}

type Hyperparams struct {
	// Diameter of transient in milliseconds.
	// Amount of time around the detected transient, for which the signal
	// will be unscaled.
	//
	// Please note that if this value is less than synthesis hop size
	// (1024 samples at 48 kHz sample rate, 21.3 ms),transient resets are
	// not guaranteed.
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
	n.detector = DetectorNew(1024*w, samplerate, opts.TransientMs, opts.PickingMs)

	return
}

func (n *Nanowarp) Process(lin, rin, lout, rout []float64, phasor *Curve) {
	fmt.Fprintln(os.Stderr, "(*Nanowarp).Process: DELETEME")

	ons := make([]float64, len(lin))
	ons1 := make([]float64, len(lin))

	if n.opts.Quality > -2 {
		poolstretch := 1.
		stretch := phasor.Dx(phasor.elems[len(phasor.elems)-1].I)
		if n.opts.ScalePool || stretch < 1 {
			poolstretch = stretch
		}
		sam := n.detector.process2(lin, rin, ons, ons1, poolstretch)

		c := phasor.Clone()
		n.bendPhasor(phasor, c, sam)
		phasor = c
	}

	n.warper.process6([][]float64{lin, rin}, [][]float64{lout, rout}, phasor)
}

func (n *Nanowarp) bendPhasor(old, new *Curve, onsets [][2]float64) {
	tsa := n.opts.TransientMs * n.fs / 1000
	for k := 0; k < len(onsets)-1; k++ {
		a := onsets[k]
		b := onsets[k+1]
		j, _ := old.Sample(b[0])
		sa := old.Dx(j)
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
		j, _ := old.Sample(i)
		a, b := new.Between(i-r), new.Between(i+r)
		new.Mutate(func(f []Breakpoint) []Breakpoint {
			if a != b {
				f = slices.Delete(f, a, b)
			}
			f = slices.Insert(f, a+1, []Breakpoint{{i - r, j - r}, {i + r, j + r}}...)
			return f
		})
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
