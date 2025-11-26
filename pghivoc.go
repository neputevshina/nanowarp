package nanowarp

import (
	"math"
	"math/rand/v2"
)

const pipi = 2 * math.Pi

func (n *Nanowarp) pghivoc(stretch float64, out []complex128) {
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
			const1 := 2.0 * math.Pi * aana / float64(n.nfft)
			const2 := 2.0 * math.Pi * asyn / float64(n.nfft)
			l := princarg(phia(w, t+1)-phia(w, t)-const1*w) / (2.0 * aana)
			r := princarg(phia(w, t)-phia(w, t-1)-const1*w) / (2.0 * aana)

			return asyn*(l-r) + const2*w
		}
		return (tau(w, t) + tau(w, t+1)) / 2
	}
	fcent := func(w, t float64) float64 {
		omega := func(w, t float64) float64 {
			a := phia(w, t) - phia(w-1, t)
			return princarg(a) / bana
		}
		if w == 0 || int(w) == n.nbins-1 {
			return phia(w, t)
		}
		return (omega(w+1, t) + omega(w, t)) / 2
	}
	for i := 0; len(n.iset) > 0 && i < n.nbins; i++ {
		h := heappop(&n.heap)
		w := float64(h.w)

		if h.t == -1 {
			if _, ok := n.iset[h.w]; ok {
				a.Phase[h.w] = a.Phase[h.w] + (tcent(w, -1)+tcent(w, 0))*asyn/2
				// a.Phase[h.w] = a.Phase[h.w] + tcent(w, 0)*asyn
				delete(n.iset, h.w)
				heappush(&n.heap, heaptriple{maga(w, 0), h.w, 0})
			}
		}
		if h.t == 0 {
			if _, ok := n.iset[h.w+1]; ok {
				a.Phase[h.w+1] = a.Phase[h.w] + (fcent(w+1, 0)+fcent(w, 0))*bsyn/2
				delete(n.iset, h.w+1)
				heappush(&n.heap, heaptriple{maga(w, 0), h.w + 1, 0})
			}
			if _, ok := n.iset[h.w-1]; ok {
				a.Phase[h.w-1] = a.Phase[h.w] - (fcent(w-1, 0)+fcent(w, 0))*bsyn/2
				delete(n.iset, h.w-1)
				heappush(&n.heap, heaptriple{maga(w, 0), h.w - 1, 0})
			}
		}

	}

	for i := range out {
		out[i] = poltocar(maga(float64(i), 0), a.Phase[i])
	}
}

// A stripped-down version of the standard container/heap.

type heaptriple struct {
	mag  float64
	w, t int
}

func less(h *[]heaptriple, i, j int) bool {
	// [...] To build a priority queue, implement the Heap interface with the
	// (negative) priority as the ordering for the Less method, so Push adds
	// items while Pop removes the highest-priority item from the queue.
	return (*h)[i].mag > (*h)[j].mag
}

func Init(h *[]heaptriple) {
	// heapify
	n := len(*h)
	for i := n/2 - 1; i >= 0; i-- {
		down(h, i, n)
	}
}

func swap(h *[]heaptriple, i, j int) {
	(*h)[i], (*h)[j] = (*h)[j], (*h)[i]
}

func pop(h *[]heaptriple) heaptriple {
	n := len(*h) - 1
	l := (*h)[n]
	*h = (*h)[:n]
	return l
}

// heappush pushes the element x onto the heap.
// The complexity is O(log n) where n = h.Len().
func heappush(h *[]heaptriple, x heaptriple) {
	*h = append(*h, x)
	up(h, len(*h)-1)
}

// heappop removes and returns the minimum element (according to Less) from the heap.
// The complexity is O(log n) where n = h.Len().
// heappop is equivalent to [heapremove](h, 0).
func heappop(h *[]heaptriple) heaptriple {
	n := len(*h) - 1
	swap(h, 0, n)
	down(h, 0, n)
	return pop(h)
}

// heapremove removes and returns the element at index i from the heap.
// The complexity is O(log n) where n = h.Len().
func heapremove(h *[]heaptriple, i int) any {
	n := len(*h) - 1
	if n != i {
		swap(h, i, n)
		if !down(h, i, n) {
			up(h, i)
		}
	}
	return pop(h)
}

// heapfix re-establishes the heap ordering after the element at index i has changed its value.
// Changing the value of the element at index i and then calling heapfix is equivalent to,
// but less expensive than, calling [heapremove](h, i) followed by a Push of the new value.
// The complexity is O(log n) where n = h.Len().
func heapfix(h *[]heaptriple, i int) {
	if !down(h, i, len(*h)) {
		up(h, i)
	}
}

func up(h *[]heaptriple, j int) {
	for {
		i := (j - 1) / 2 // parent
		if i == j || less(h, j, i) {
			break
		}
		swap(h, i, j)
		j = i
	}
}

func down(h *[]heaptriple, i0, n int) bool {
	i := i0
	for {
		j1 := 2*i + 1
		if j1 >= n || j1 < 0 { // j1 < 0 after int overflow
			break
		}
		j := j1 // left child
		if j2 := j1 + 1; j2 < n && less(h, j2, j1) {
			j = j2 // = 2*i + 2  // right child
		}
		if !less(h, j, i) {
			break
		}
		swap(h, i, j)
		i = j
	}
	return i > i0
}
