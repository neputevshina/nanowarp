package nanowarp

import (
	"fmt"
	"slices"
)

// ErrNonMonotonicCurve is returned when the given [Curve] is not monotonic when a monotonic one is expected.
type ErrNonMonotonicCurve struct {
	Index int
}

func (i *ErrNonMonotonicCurve) Error() string {
	return fmt.Sprintf(`curve is not monotonic, e[%d+1].J<e[%d].J`, i.Index, i.Index)
}

// Onset is an onset point of a transient.
type Onset struct {
	I     float64 // Input sample index.
	Power float64 // Absolute unitless power of a transient.
}

// Breakpoint is a point of a [Curve].
//
// I is the input sample index, and J is the output sample index.
type Breakpoint struct {
	I, J float64
}

// Bp is a quick constructor for a [Breakpoint].
func Bp(i, j float64) Breakpoint { return Breakpoint{I: i, J: j} }

// Curve is a line-segment curve, describing the time mapping between input sample indices
// and output sample indices.
//
// Output sample indices are used as an input to curve function and Curve is always guaranteed to
// be indexable by J.
//
// A Curve is indexable by I iff it is monotonic.
// Currently, Nanowarp does not support variable stretching using non-monotonic curves.
type Curve struct {
	elems       []Breakpoint
	last, rlast int
	start, end  Breakpoint
}

// NewCurve creates a new [Curve] from individual breakpoints.
//
// The curve must be monotonic.
// If not, [ErrNonMonotonicCurve] is returned.
func NewCurve(bps []Breakpoint) (*Curve, error) {
	c := &Curve{}
	c.Mutate(func(b []Breakpoint) []Breakpoint {
		return slices.Clone(bps)
	})
	return c, c.Validate()
}

// Dx returns the derivative of the curve with respect to input sample scale
// at the given input sample offset.
func (c *Curve) Dx(i float64) (v float64) {
	if i >= c.end.I {
		return c.dx(len(c.elems) - 2)
	}
	if i < c.start.I {
		return c.dx(0)
	}
	f := c.Between(i)
	return c.dx(f)
}

// Dx returns the derivative of the curve with respect to output sample scale
// at the given output sample offset.
func (c *Curve) Dy(j float64) (v float64) {
	if j >= c.end.J {
		return 1 / c.dx(len(c.elems)-2)
	}
	if j < c.start.J {
		return 1 / c.dx(0)
	}
	f := c.ReverseBetween(j)
	return 1 / c.dx(f)
}

func (c *Curve) dx(f int) float64 {
	delx := (c.elems[f+1].I - c.elems[f].I)
	dely := (c.elems[f+1].J - c.elems[f].J)
	return dely / delx
}

// Sample returns the value of the curve at the given input sample index.
// This function is not guaranteed to return a correct result if the curve
// is not monotonic.
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

// ReverseSample returns the value of the curve at the given output sample index.
func (c *Curve) ReverseSample(j float64) (i float64) {
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
func (c *Curve) Between(i float64) (a int) {
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
func (c *Curve) ReverseBetween(j float64) (a int) {
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
// maintaining the validity of a Curve object.
func (c *Curve) Mutate(f func([]Breakpoint) []Breakpoint) {
	c.elems = f(c.elems)
	c.mutate()
}

func (c *Curve) mutate() {
	c.start = c.elems[0]
	c.end = c.elems[len(c.elems)-1]
	c.last = 0
	c.rlast = 0
}

// Clone returns the copy of a Curve.
func (c *Curve) Clone() *Curve {
	oc := &Curve{
		elems: slices.Clone(c.elems),
	}
	oc.mutate()
	return oc
}

// Validate checks the curve for monotonicity and if it is
// not monotonic returns [ErrNonMonotonicCurve].
// It returns nil otherwise.
func (c *Curve) Validate() error {
	for e := range c.elems[1:] {
		if c.elems[e+1].J < c.elems[e].J {
			return &ErrNonMonotonicCurve{Index: e}
		}
	}
	return nil
}
