package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
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
	for i := e; i < len(in); i += e {
		enfft := func(i int, x []complex128, w []float64) {
			zero(a.S)
			copy(a.S, in[i:i+n.nbuf])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)
		enfft(i, a.Xd, a.Wd)
		enfft(i, a.Xt, a.Wt)

		if i == e {
			copy(a.P, a.X)
			continue
		}

		// Begin of PGHI

		a := &n.a
		n.heap = n.heap[:0]
		clear(n.iset)

		abstol := mag(a.X[0])
		for j := range a.X {
			abstol = max(abstol, mag(a.X[j]), mag(a.P[j]))
		}
		abstol *= 1e-6

		for j := range a.X {
			if mag(a.X[j]) > abstol {
				n.iset[j] = struct{}{}
				heappush(&n.heap, heaptriple{mag(a.X[j]), j, -1})
			} else {
				a.Phase[j] = mix(-math.Pi, math.Pi, rand.Float64())
			}
		}

		for j := range a.X {
			var f float64 = float64(j) / float64(len(a.X))
			_ = f
			// a.S2[j] = princarg(f + math.Pi + -2*imag(a.Xd[j]/a.X[j]))
			// a.S2[j] = -imag(a.Xd[j]/a.X[j]) * (0.001 * float64(i/e))
			// a.S2[j] = -imag(a.Xd[j]/a.X[j]) * (0.001*16000 + 0.001*float64(i/e))
			// a.S2[j] = princarg(f*math.Pi - imag(a.Xd[j]/a.X[j])*(0.001*float64(i/e)))
			// a.S2[j] = princarg(f*math.Pi - imag(a.Xd[j]/a.X[j])*(8+0.001*float64(i/e)))
			// a.S2[j] = princarg(2*f - 2*imag(a.Xd[j]/a.X[j]))
			// a.S2[j] = princarg(0.01*float64(i/e) - imag(a.Xd[j]/a.X[j]))
			a.S2[j] = princarg(0.75*float64(j) + 0.25*imag(a.Xd[j]/a.X[j]))
			// a.S2[j] = cmplx.Phase(a.X[j])
			// a.S2[j] = princarg(cmplx.Phase(a.X[j]) - cmplx.Phase(a.Cs0[j]))
		}
		G[`phasogram.png`] = append(G[`phasogram.png`].([][]float64), slices.Clone(a.S2))

		for j := range a.X {
			var f float64 = float64(j) / float64(len(a.X))
			_ = f
			a.Phase[j] = princarg(cmplx.Phase(a.P[j]) + a.S2[j])
			// a.Phase[j] = mix(-math.Pi, math.Pi, rand.Float64())
		}

		// for i := 0; len(n.iset) > 0 && i < n.nbins; i++ {
		// 	h := heappop(&n.heap)
		// 	dt := a.S2[h.w] * stretch
		// 	dw := (math.Pi - real(a.Xt[i]/a.X[i])/float64(len(a.X))*2*math.Pi) * stretch
		// 	switch h.t {
		// 	case -1:
		// 		if _, ok := n.iset[h.w]; ok {
		// 			a.Phase[h.w] = princarg(cmplx.Phase(a.P[h.w]) + dt)
		// 			delete(n.iset, h.w)
		// 			heappush(&n.heap, heaptriple{mag(izero(a.X, h.w)), h.w, 0})
		// 		}

		// 	case 0:
		// 		if _, ok := n.iset[h.w+1]; ok {
		// 			a.Phase[h.w+1] = princarg(a.Phase[h.w] + dw)
		// 			delete(n.iset, h.w+1)
		// 			heappush(&n.heap, heaptriple{mag(izero(a.X, h.w)), h.w + 1, 0})
		// 		}
		// 		if _, ok := n.iset[h.w-1]; ok {
		// 			a.Phase[h.w-1] = princarg(a.Phase[h.w] - dw)
		// 			delete(n.iset, h.w-1)
		// 			heappush(&n.heap, heaptriple{mag(izero(a.X, h.w)), h.w - 1, 0})
		// 		}
		// 	}
		// }

		// G[`phasogram.png`] = append(G[`phasogram.png`].([][]float64), slices.Clone(a.Phase))

		for j := range a.Phase {
			a.A[j] = cmplx.Rect(mag(a.X[j]), a.Phase[j])
		}

		copy(a.P, a.A)
		copy(a.Cs0, a.X)

		// End of PGHI

		n.fft.Sequence(a.S, a.P)
		for j := range a.S {
			a.S[j] /= n.norm
			// a.S[i] *= stretch
		}
		mul(a.S, a.W)
		add(out[adv:adv+n.nbuf], a.S)
		adv += o
	}
}
