package nanowarp

import (
	"math"
	"math/rand/v2"
)

func (n *Nanowarp) pghi(stretch float64, out []complex128) {
	a := &n.a
	n.heap = n.heap[:0]
	clear(n.iset)

	abstol := mag(a.X[0])
	for i := range a.X {
		abstol = max(abstol, mag(a.X[i]), mag(a.P[i]))
	}
	abstol *= 1e-6

	G[`phasogram.png`] = append(G[`phasogram.png`].([][]float64), make([]float64, len(a.Phase)))
	G[`origphase.png`] = append(G[`origphase.png`].([][]float64), make([]float64, len(a.Phase)))
	G[`mag.png`] = append(G[`mag.png`].([][]float64), make([]float64, len(a.Phase)))

	for i := range a.X {
		if mag(a.X[i]) > abstol {
			n.iset[i] = struct{}{}
			heappush(&n.heap, heaptriple{mag(a.P[i]), i, -1})
		} else {
			a.Phase[i] = mix(-math.Pi, math.Pi, rand.Float64())
		}
	}

	aana := float64(n.hop)
	asyn := aana * stretch
	bana := float64(n.hop)
	bsyn := bana * stretch
	for i := 0; len(n.iset) > 0 && i < n.nbins; i++ {
		h := heappop(&n.heap)
		w := float64(h.w)
		phia := func(w, t float64) float64 {
			t += 2
			return angle(izero(a.F[int(t)], int(w)))
		}
		tcent := func(w, t float64) float64 {
			tau := func(w, t float64) float64 {
				ramp := pipi * w / float64(n.nfft)
				a := phia(w, t) - phia(w, t-1) - aana*ramp
				return princarg(a)/aana + ramp
			}
			return (tau(w, t) + tau(w, t+1)) / 2
		}
		fcent := func(w, t float64) float64 {
			omega := func(w, t float64) float64 {
				a := phia(w, t) - phia(w-1, t)
				return princarg(a) / bana
			}
			return (omega(w+1, t) + omega(w, t)) / 2
		}

		if h.t == -1 {
			if _, ok := n.iset[h.w]; ok {
				a.Phase[h.w] = a.Phase[h.w] + (tcent(w, -1)+tcent(w, 0))*asyn/2
				delete(n.iset, h.w)
				heappush(&n.heap, heaptriple{mag(a.X[h.w]), h.w, 0})
			}
		}
		if h.t == 0 {
			if _, ok := n.iset[h.w+1]; ok {
				a.Phase[h.w+1] = a.Phase[h.w] + (fcent(w, -1)+fcent(w, 0))*bsyn/2
				delete(n.iset, h.w+1)
				heappush(&n.heap, heaptriple{mag(a.X[h.w+1]), h.w + 1, 0})
			}
			if _, ok := n.iset[h.w-1]; ok {
				a.Phase[h.w-1] = a.Phase[h.w] - (fcent(w, -1)+fcent(w, 0))*bsyn/2
				delete(n.iset, h.w-1)
				heappush(&n.heap, heaptriple{mag(a.X[h.w-1]), h.w - 1, 0})
			}
		}

	}

	for i := range out {
		out[i] = poltocar(mag(a.X[i]), a.Phase[i])
	}
}
