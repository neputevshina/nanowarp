package nanowarp

import (
	"math"
	"slices"
	"sync"

	"github.com/neputevshina/nanowarp/dspio"
)

// Nanowarp is a time-scale modification algorithm.
type Nanowarp struct {
	fs, nch int

	warper   *warper
	detector *detector

	opts      Options
	stretch   float64
	semitones float64
}

// Options contains configuration options for Nanowarp.
type Options struct {
	// Output scaled onsets only.
	Onsets bool

	// Set algorithm quality.
	//  -1: Use brute force approximation to PGHI. Less transparent, 20% faster.
	//  0:  Use PGHI.
	Quality int

	// Enable phase resets.
	//  -2: Don't perform transient separation, output raw PVDR without phase resets.
	//  -1: Use original phase when not stretching.
	//      Introduces clicky artifacts but cleanest for transient-heavy material.
	//      Best numerical stability because of resets.
	//  0:  When not stretching, advance phase for harmonic components,
	//      use original phase otherwise.
	//      Very little to no artifacts, but noticeable slight loss in clarity.
	Resets int

	// Channel for receiving processing progress.
	//
	// If not nil, this channel will receive current input and output sample index
	// pair for every 5 seconds of output and at the start and end of processing.
	//
	// Nanowarp will close the channel after the end of processing.
	Progress chan<- Progress

	Hyperparams
}

// Progress is a message containing the progress of an algorithm.
type Progress struct {
	Current float64
	End     float64
	Process string
}

// Hyperparams contains low-level configuration of Nanowarp algorithm.
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
	// Limits vertical propagation of ridges' region of influence.
	//
	// Region of influence a ridge is defined as a set of all linked arrow
	// chains starting from it, that are not ridges itself, and limited by a sink
	// (bin that have no arrows exiting from it).
	// Phase never be reset at this number of bins around the ridge.
	//
	// Higher values compromise transient quality over tonal quality.
	InfluenceRadius int `default:"3"`

	// Disable limiting of local group delay.
	//
	// Disables fix for ridge triplication in time, which is obvious on
	// extreme (>4x) stretches. Makes no effect for shrinks.
	//
	// On small coefficients result may sound more “full” and less “plastic”.
	NoTriplicationFix bool
}

// New creates a new time-scale modification process object.
func New(samplerate int, channels int, opts Options) (n *Nanowarp) {
	structinit(&opts)
	structinit(&opts.Hyperparams)
	n = new(channels, samplerate, &opts)
	n.opts = opts

	return
}

func new(channels int, samplerate int, opts *Options) (n *Nanowarp) {
	// TODO Fixed absolute bandwidth through zero-padding.
	// Hint: nbuf is already there.
	// TODO Find optimal bandwidths.
	n = &Nanowarp{}
	n.fs = samplerate
	n.nch = channels
	scale := func(x float64) int {
		e := int(math.Ceil(x * float64(samplerate) / 48000))
		return e - e%warperOverlap + warperOverlap // Quantize by warper overlap.
	}

	n.warper = warperNew(scale(4096), 2, channels, n)
	n.detector = detectorNew(scale(1024), samplerate, channels, opts.PickingMs)

	return
}

// Lengther returns the size of the underlying stream in amount
// of multichannel samples.
type Lengther interface {
	Length() int
}

// Process pefrorms the time-scale modification of r and writes it to w.
func (n *Nanowarp) Process(r dspio.GrainReadSeeker, w dspio.SignalWriter, phasor *Curve) {
	if n.opts.Resets > -2 {
		poolstretch := 1.
		stretch := phasor.Dx(phasor.elems[len(phasor.elems)-1].I)
		if n.opts.ScalePool || stretch < 1 {
			poolstretch = stretch
		}
		po, pi := dspio.GoPipe(2)
		wg := sync.WaitGroup{}
		wg.Add(3)
		sam := make([]Onset, 0)
		onsc := make(chan Onset, 0)
		go func() {
			defer wg.Done()
			err := n.detector.noveltyCurveProcess(r, pi)
			if err != nil {
				panic(err)
			}
			pi.Close()
		}()
		go func() {
			defer wg.Done()
			err := n.detector.dilatePeakSelectProcess(po, nil, poolstretch, onsc)
			if err != nil {
				panic(err)
			}
		}()
		go func() {
			filelen := 0.

			l, ok := r.(Lengther)
			if ok {
				filelen = float64(l.Length())
			}

			defer wg.Done()
			for o := range onsc {
				if n.opts.Progress != nil {
					if !ok {
						filelen = o.I
					}
					n.opts.Progress <- Progress{
						Current: o.I,
						End:     filelen,
						Process: "Analysis",
					}
				}
				sam = append(sam, o)
			}
		}()
		wg.Wait()

		c := phasor.Clone()
		n.bendPhasor(phasor, c, sam)

		phasor = c
	}

	grw := dspio.NewRegularToOfflineGrainWriter(n.warper.nbuf, n.warper.hop, w)
	n.warper.process(r, grw, phasor)
}

func (n *Nanowarp) bendPhasor(old, new *Curve, onsets []Onset) {
	tsa := n.opts.TransientMs * n.fs / 1000
	for k := 0; k < len(onsets)-1; k++ {
		a := onsets[k]
		b := onsets[k+1]
		j, _ := old.Sample(b.I)
		sa := old.Dx(j)
		if b.I-a.I < float64(max(n.warper.nbuf, tsa))/sa {
			// Leave only louder one
			if a.Power > b.Power {
				onsets[k] = a
			}
			copy(onsets[k:], onsets[k+1:])
			onsets = onsets[:len(onsets)-1]
			k--
		}
	}
	for k := 0; k < len(onsets)-1; k++ {
		i := onsets[k].I
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
