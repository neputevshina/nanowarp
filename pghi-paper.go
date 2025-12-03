package nanowarp

import (
	"container/heap"
	"math"
	"math/rand/v2"
	"slices"
)

func (n *Nanowarp) pghipaper(stretch float64, out []complex128) {
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
			heap.Push(&n.heap, heaptriple{mag(a.P[i]), i, -1})
		} else {
			a.Phase[i] = mix(-math.Pi, math.Pi, rand.Float64())
		}
	}

	// aana := float64(n.hop)
	// asyn := aana * stretch
	// bana := float64(n.hop)
	// bsyn := bana * stretch
	aana := float64(n.hop) / stretch
	asyn := float64(n.hop)
	bana := float64(n.hop) / stretch
	bsyn := float64(n.hop)
	phia := func(w, t float64) float64 {
		t += 2
		return angle(izero(a.F[int(t)], int(w)))
	}
	maga := func(w, t float64) float64 {
		t += 2
		return mag(izero(a.F[int(t)], int(w)))
	}
	tcent := func(w, t float64) float64 {
		tau := func(w, t float64) float64 {
			ramp := 2 * math.Pi * w / float64(n.nfft)
			// ramp = 0
			a := phia(w, t) - phia(w, t-1) - aana*ramp
			return princarg(a)/aana + ramp
		}
		return (tau(w, t) + tau(w, t+1)) / 2
	}
	fcent := func(w, t float64) float64 {
		omega := func(w, t float64) float64 {
			a := phia(w, t) - phia(w-1, t)
			return princarg(a) / bana
			// return a / bana
		}
		if w == 0 || int(w) == n.nbins-1 {
			return phia(w, t)
		}
		return (omega(w+1, t) + omega(w, t)) / 2
	}
	for i := range a.X {
		w := float64(i)
		a.S2[i] = min(50, max(-10, (tcent(w, -1)+tcent(w, 0))*asyn/2))
		// a.S2[i] = min(4, max(-4, ((fcent(w+1, 0)+fcent(w, 0))*bsyn/2+math.Pi)))
		// a.S2[i] = mag(a.X[i])
	}
	G[`phasogram.png`] = append(G[`phasogram.png`].([][]float64), slices.Clone(a.S2))

	for i := 0; len(n.iset) > 0 && i < n.nbins; i++ {
		h := heap.Pop(&n.heap).(heaptriple)
		w := float64(h.w)

		if h.t == -1 {
			if _, ok := n.iset[h.w]; ok {
				a.Phase[h.w] = a.Phase[h.w] + (tcent(w, -1)+tcent(w, 0))*asyn/2
				// a.Phase[h.w] = a.Phase[h.w] + tcent(w, 0)*asyn
				delete(n.iset, h.w)
				heap.Push(&n.heap, heaptriple{maga(w, 0), h.w, 0})
			}
		}
		if h.t == 0 {
			if _, ok := n.iset[h.w+1]; ok {
				a.Phase[h.w+1] = a.Phase[h.w] + (fcent(w+1, 0)+fcent(w, 0))*bsyn/2
				delete(n.iset, h.w+1)
				heap.Push(&n.heap, heaptriple{maga(w, 0), h.w + 1, 0})
			}
			if _, ok := n.iset[h.w-1]; ok {
				a.Phase[h.w-1] = a.Phase[h.w] - (fcent(w-1, 0)+fcent(w, 0))*bsyn/2
				delete(n.iset, h.w-1)
				heap.Push(&n.heap, heaptriple{maga(w, 0), h.w - 1, 0})
			}
		}

	}

	for i := range out {
		out[i] = poltocar(maga(float64(i), 0), a.Phase[i])
	}
}
