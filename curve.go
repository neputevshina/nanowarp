package nanowarp

import (
	"fmt"
	"slices"
)

type ErrorNotMonotonicCurve struct {
	Index int
}

func (i *ErrorNotMonotonicCurve) Error() string {
	return fmt.Sprintf(`curve is not monotonic, e[%d+1]<e[%d]`, i.Index, i.Index)
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

func NewCurve(bps []Breakpoint) (*Curve, error) {
	c := &Curve{}
	c.Mutate(func(b []Breakpoint) []Breakpoint {
		return slices.Clone(bps)
	})
	return c, c.Validate()
}

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

func (c *Curve) IntDy(j int) (v float64) {
	return c.Dy(float64(j))
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

func (c *Curve) ReverseSample(j float64) (i float64) {
	if j >= c.end.J {
		return c.end.I
	}
	if j < c.start.J {
		return c.start.I
	}
	f := c.ReverseBetween(j)
	nj := unmix(c.elems[f].J, c.elems[f+1].J, j)
	i = precisionmix(c.elems[f].I, c.elems[f+1].I, nj)
	return
}

func (c *Curve) IntReverseSample(j int) (i int) {
	return int(c.ReverseSample(float64(j)))
}

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

func (c *Curve) Clone() *Curve {
	oc := &Curve{
		elems: slices.Clone(c.elems),
	}
	oc.mutate()
	return oc
}

func (c *Curve) Validate() error {
	for e := range c.elems[1:] {
		if c.elems[e+1].I < c.elems[e].I ||
			c.elems[e+1].J < c.elems[e].J {
			return &ErrorNotMonotonicCurve{Index: e}
		}
	}
	return nil
}
