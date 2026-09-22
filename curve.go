package nanowarp

import (
	"fmt"
	"slices"
)

// Tween describes the curve type between two [Phasorpoint]s,
// with former taking the lead.
type Tween int

// Tween types.
const (
	// Zero-order hold, same value for the whole range.
	// For [Phasorpoint] is the same as [TweenLinear].
	TweenZoh Tween = iota
	// Linear segment.
	TweenLinear
)

// ErrNonMonotonicPhasor is the error returned when the given [Phasor] is not monotonic when a monotonic one is expected.
type ErrNonMonotonicPhasor struct {
	Index int
}

// Error implements [error].
func (i *ErrNonMonotonicPhasor) Error() string {
	return fmt.Sprintf(`curve is not monotonic, e[%d+1].J<e[%d].J`, i.Index, i.Index)
}

// ErrInvalidCurve is the error returned when the given [Phasor] or [Envelope] has a point when J is less than J of a previous point.
type ErrInvalidCurve struct {
	Index int
}

// Error implements [error].
func (i *ErrInvalidCurve) Error() string {
	return fmt.Sprintf(`curve is not monotonic, e[%d+1].J<e[%d].J`, i.Index, i.Index)
}

// Onset is an onset point of a transient.
type Onset struct {
	I     float64 // Input sample index.
	Power float64 // Absolute unitless power of a transient.
}

// Phasorpoint is a point of a [Phasor].
//
// J is the output sample index and I is the input sample index.
type Phasorpoint struct {
	J, I float64
	Tween

	reset bool
}

// Pp is a quick constructor for a [Phasorpoint].
func Pp(j, i float64) Phasorpoint { return Phasorpoint{I: i, J: j} }

// Breakpoint is a point of an [Envelope].
//
// J is the sample index and V is the value.
type Breakpoint struct {
	J, V float64
	Tween
}

// Bp is a quick constructor for a [Breakpoint].
func Bp(j, value float64) Breakpoint { return Breakpoint{J: j, V: value} }

// Phasor is a curve describing the time mapping between input sample indices and output sample indices.
//
// Output sample indices are used as an input to curve function and a valid Phasor is always guaranteed to
// be indexable by J.
//
// A Phasor is indexable by I iff it is monotonic.
//
// Currently, Nanowarp does not support variable stretching using non-monotonic curves (can't reverse while stretching).
type Phasor struct {
	elems       []Phasorpoint
	last, rlast int
	start, end  Phasorpoint
}

// NewPhasor creates a new [Phasor] from individual breakpoints.
//
// The curve must be monotonic.
// If not, [ErrNonMonotonicPhasor] is returned.
func NewPhasor(bps []Phasorpoint) (*Phasor, error) {
	c := &Phasor{}
	c.Mutate(func(b []Phasorpoint) []Phasorpoint {
		return slices.Clone(bps)
	})
	return c, c.Validate()
}

// Dx returns the derivative of the curve with respect to input sample scale
// at the given input sample offset.
func (c *Phasor) Dx(i float64) (v float64) {
	if i >= c.end.I {
		return c.dx(len(c.elems) - 2)
	}
	if i < c.start.I {
		return c.dx(0)
	}
	f := c.Between(i)
	return c.dx(f)
}

// Dy returns the derivative of the curve with respect to output sample scale
// at the given output sample offset.
func (c *Phasor) Dy(j float64) (v float64) {
	if j >= c.end.J {
		return 1 / c.dx(len(c.elems)-2)
	}
	if j < c.start.J {
		return 1 / c.dx(0)
	}
	f := c.ReverseBetween(j)
	return 1 / c.dx(f)
}

func (c *Phasor) dx(f int) float64 {
	delx := (c.elems[f+1].I - c.elems[f].I)
	dely := (c.elems[f+1].J - c.elems[f].J)
	return dely / delx
}

// Sample returns the value of a curve at the given input sample index.
// This function is not guaranteed to return a correct result if the curve
// is not monotonic.
func (c *Phasor) Sample(i float64) (j float64, oflow int) {
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

// ReverseSample returns the value of the curve at the given output sample index.
func (c *Phasor) ReverseSample(j float64) (i float64) {
	if j >= c.end.J {
		return c.end.I + j - c.end.J
	}
	if j < c.start.J {
		return c.start.I + j - c.start.J
	}
	f := c.ReverseBetween(j)
	nj := unmix(c.elems[f].J, c.elems[f+1].J, j)
	i = precisionmix(c.elems[f].I, c.elems[f+1].I, nj)
	return
}

// Between returns an integer index of an internal slice of Breakpoints such as
// sl[a].I < i < sl[a+1].I.
func (c *Phasor) Between(i float64) (a int) {
	if c.elems[c.last].I < i {
		c.last = 0
	}
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

// ReverseBetween returns an integer index of an internal slice of Breakpoints such as
// sl[a].J < j < sl[a+1].J.
func (c *Phasor) ReverseBetween(j float64) (a int) {
	if c.elems[c.rlast].I < j {
		c.rlast = 0
	}
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

// Mutate allows editing the internal slice of Breakpoints,
// maintaining the validity of a Phasor object.
func (c *Phasor) Mutate(f func([]Phasorpoint) []Phasorpoint) {
	c.elems = f(c.elems)
	c.mutate()
}

func (c *Phasor) mutate() {
	c.start = c.elems[0]
	c.end = c.elems[len(c.elems)-1]
	c.last = 0
	c.rlast = 0
}

// Clone returns the copy of a Curve.
func (c *Phasor) Clone() *Phasor {
	oc := &Phasor{
		elems: slices.Clone(c.elems),
	}
	oc.mutate()
	return oc
}

// Validate checks the curve for correctness and returns [ErrInvalidCurve] when
// it is invalid or [ErrNonMonotonicPhasor] when it is not monotonic.
// It returns nil otherwise.
func (c *Phasor) Validate() error {
	for e := range c.elems[1:] {
		if c.elems[e+1].J < c.elems[e].J {
			return &ErrInvalidCurve{Index: e}
		}
		if c.elems[e+1].I < c.elems[e].I {
			return &ErrNonMonotonicPhasor{Index: e}
		}
	}
	return nil
}

// reverseReset returns the state of a reset flag for a given j.
func (c *Phasor) reverseReset(j float64) bool {
	if j >= c.end.J || j < c.start.J {
		return false
	}
	f := c.ReverseBetween(j)
	return c.elems[f].reset
}

// Envelope is a curve describing the value of a variable mapped by sample indices.
type Envelope struct {
	elems       []Breakpoint
	last, rlast int
	start, end  Breakpoint
}

// NewEnvelope creates a new [Envelope] from individual breakpoints.
func NewEnvelope(bps []Breakpoint) (*Envelope, error) {
	c := &Envelope{}
	c.Mutate(func(b []Breakpoint) []Breakpoint {
		return slices.Clone(bps)
	})
	return c, c.Validate()
}

// Dy returns the derivative of the curve value at the given sample offset.
func (c *Envelope) Dy(j float64) (v float64) {
	if j >= c.end.J {
		return c.dx(len(c.elems) - 2)
	}
	if j < c.start.J {
		return c.dx(0)
	}
	f := c.Between(j)
	return c.dx(f)
}

func (c *Envelope) dx(f int) float64 {
	dela := (c.elems[f+1].V - c.elems[f].V)
	dely := (c.elems[f+1].J - c.elems[f].J)
	return dely / dela
}

// Sample returns the value of the curve at the given output sample index.
func (c *Envelope) Sample(j float64) (v float64) {
	if j >= c.end.J {
		return c.end.V
	}
	if j < c.start.J {
		return c.start.V
	}
	f := c.Between(j)
	nj := unmix(c.elems[f].J, c.elems[f+1].J, j)
	switch c.elems[f].Tween {
	case TweenZoh:
		v = c.elems[f].V
	case TweenLinear:
		v = precisionmix(c.elems[f].V, c.elems[f+1].V, nj)
	default:
		panic(`invalid Tween type`)
	}
	return
}

// Between returns an integer index of an internal slice of Breakpoints such as
// sl[a].J < j < sl[a+1].J.
func (c *Envelope) Between(j float64) (a int) {
	if c.elems[c.rlast].J < j {
		c.rlast = 0
	}
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

// Mutate allows editing the internal slice of Breakpoints,
// maintaining the validity of an Envelope object.
func (c *Envelope) Mutate(f func([]Breakpoint) []Breakpoint) {
	c.elems = f(c.elems)
	c.mutate()
}

func (c *Envelope) mutate() {
	c.start = c.elems[0]
	c.end = c.elems[len(c.elems)-1]
	c.last = 0
	c.rlast = 0
}

// Clone returns the copy of a Curve.
func (c *Envelope) Clone() *Envelope {
	oc := &Envelope{
		elems: slices.Clone(c.elems),
	}
	oc.mutate()
	return oc
}

// Validate checks the curve for correctness and returns [ErrInvalidCurve] when
// it is invalid.
// It returns nil otherwise.
func (c *Envelope) Validate() error {
	for e := range c.elems[1:] {
		if c.elems[e+1].J < c.elems[e].J {
			return &ErrInvalidCurve{Index: e}
		}
	}
	return nil
}
