package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"slices"
)

func (n *Nanowarp) Process1(in []float64, out []float64, stretch float64) {
	a := &n.a
	e := int(math.Floor(float64(n.hop)))
	o := int(math.Floor(float64(n.hop) * stretch))
	fmt.Fprintln(os.Stderr, `outhop:`, o)
	fmt.Fprintln(os.Stderr, `inhop:`, e, `nbuf:`, n.nbuf)

	adv := 0
	for i := 0; i < len(in); i += e {
		enfft := func(i int, x []complex128, w []float64) {
			zero(a.S)
			copy(a.S, in[i:i+n.nbuf])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)
		enfft(i, a.Xd, a.Wd)
		enfft(i, a.Xt, a.Wt)

		// Begin of PGHI

		a := &n.a
		n.heap = n.heap[:0]
		clear(n.iset)

		abstol := mag(a.X[0])
		for i := range a.X {
			abstol = max(abstol, mag(a.X[i]), mag(a.P[i]))
		}
		abstol *= 1e-6

		for i := range a.X {
			if mag(a.X[i]) > abstol {
				n.iset[i] = struct{}{}
				heappush(&n.heap, heaptriple{mag(a.P[i]), i, -1})
			} else {
				// a.Phase[i] = mix(-math.Pi, math.Pi, rand.Float64())
			}
		}

		for i := range a.X {
			a.S2[i] = princarg(real(a.Xt[i]/a.X[i]) / float64(len(a.X)) * 2 * math.Pi)
			// a.S2[i] = min(4, max(-4, -imag(a.Xd[i]/a.X[i])))
			// a.S2[i] = mag(a.X[i])
		}
		G[`phasogram.png`] = append(G[`phasogram.png`].([][]float64), slices.Clone(a.S2))

		for i := 0; len(n.iset) > 0 && i < n.nbins; i++ {
			h := heappop(&n.heap)
			dt := -imag(a.Xd[h.w]/a.X[h.w]) * stretch
			dw := real(a.Xt[i]/a.X[i]) / float64(len(a.X)) * 2 * math.Pi * stretch
			switch h.t {
			case -1:
				if _, ok := n.iset[h.w]; ok {
					a.Phase[h.w] = princarg(cmplx.Phase(a.P[h.w]) + dt)
					delete(n.iset, h.w)
					heappush(&n.heap, heaptriple{mag(izero(a.X, h.w)), h.w, 0})
				}

			case 0:
				if _, ok := n.iset[h.w+1]; ok {
					a.Phase[h.w+1] = princarg(a.Phase[h.w] + dw)
					delete(n.iset, h.w+1)
					heappush(&n.heap, heaptriple{mag(izero(a.X, h.w)), h.w + 1, 0})
				}
				if _, ok := n.iset[h.w-1]; ok {
					a.Phase[h.w-1] = princarg(a.Phase[h.w] - dw)
					delete(n.iset, h.w-1)
					heappush(&n.heap, heaptriple{mag(izero(a.X, h.w)), h.w - 1, 0})
				}
			}
		}

		for i := range a.Phase {
			a.A[i] = cmplx.Rect(mag(a.X[i]), a.Phase[i])
		}

		copy(a.P, a.A)

		// End of PGHI

		n.fft.Sequence(a.S, a.A)
		for i := range a.S {
			a.S[i] /= n.norm
			// a.S[i] *= stretch
		}
		mul(a.S, a.W)
		add(out[adv:adv+n.nbuf], a.S)
		adv += o
	}
}
