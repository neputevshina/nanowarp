package nanowarp

import (
	"math"
	"math/rand/v2"
)

const pipi = 2 * math.Pi

func (n *Nanowarp) pghi(stretch float64, out []complex128) {
	// TODO:
	// - Remove phase math, replace it with complex.
	//	- For frequency reassignment use sin and cos. Will be faster than atan2.
	// - Pre-compute mag(x) and mag(xp).
	// - Pre-compute random phase texture.
	a := &n.a
	n.heap = n.heap[:0]
	clear(n.iset)
	// aap := int(math.Floor(float64(n.hop) * n.prevstretch))
	// aac := int(math.Floor(float64(n.hop) * stretch))

	abstol := mag(a.X[0])
	for i := range a.X {
		abstol = max(abstol, mag(a.X[i]), mag(a.P[i]))
	}
	abstol *= 1e-6

	// for i := range a.X {
	// 	a.Frame[i] = 0
	// }

	G[`phasogram`] = append(G[`phasogram`].([][]float64), make([]float64, len(a.Frame)))
	G[`origphase`] = append(G[`origphase`].([][]float64), make([]float64, len(a.Frame)))
	G[`deltaphase`] = append(G[`deltaphase`].([][]float64), make([]float64, len(a.Frame)))
	G[`sigmaphase`] = append(G[`sigmaphase`].([][]float64), make([]float64, len(a.Frame)))

	for i := range a.X {
		if mag(a.X[i]) > abstol {
			n.iset[i] = struct{}{}
			heappush(&n.heap, heaptriple{mag(a.P[i]), i, -1})
		} else {
			// FIXME Slow. Use noise texture.
			a.Frame[i] = mix(-math.Pi, math.Pi, rand.Float64())
		}
	}

	aana := float64(n.hop)
	asyn := aana * stretch
	for i := range a.X {
		if i == n.nbins-1 || i == 0 {
			a.S2[i] = angle(a.X[i])
		} else {
			a.S2[i] = (princarg(angle(a.X[i+1])-angle(a.X[i])) +
				princarg(angle(a.X[i])-angle(a.X[i-1]))) / 2
		}
		const1 := 2.0 * math.Pi * aana / float64(n.nfft)
		const2 := 2.0 * math.Pi * asyn / float64(n.nfft)
		f := func(j int) float64 {
			return asyn*(princarg(angle(a.F[j+1][i])-angle(a.F[j][i])-const1*float64(i))/(2.0*aana)-
				princarg(angle(a.F[j][i])-angle(a.F[j-1][i])-const1*float64(i))/(2.0*aana)) + const2*float64(i)
		}
		a.S3[i] = f(1)
		a.S4[i] = f(2)
		// 	return -real(a.Xt[i] / (a.X[i] + eps)) / float64(nbins) * 2 * math.Pi

	}
	for i := 0; len(n.iset) > 0 && i < n.nbins; i++ {
		h := heappop(&n.heap)
		dt := func(j int) float64 {
			if j == -1 {
				return a.S3[h.w]
			}
			return a.S4[h.w]
		}
		df := func(i int) float64 {
			return a.S2[i]
		}

		if h.t == -1 {
			if _, ok := n.iset[h.w]; ok {
				a.Frame[h.w] = a.Frame[h.w] + (dt(-1)+dt(0))/2
				delete(n.iset, h.w)
				heappush(&n.heap, heaptriple{mag(a.X[h.w]), h.w, 0})
			}
		}
		if h.t == 0 {
			if _, ok := n.iset[h.w+1]; ok {
				a.Frame[h.w+1] = a.Frame[h.w] + (df(h.w)+df(h.w+1))/2*stretch
				delete(n.iset, h.w+1)
				heappush(&n.heap, heaptriple{mag(a.X[h.w+1]), h.w + 1, 0})
			}
			if _, ok := n.iset[h.w-1]; ok {
				a.Frame[h.w-1] = a.Frame[h.w] - (df(h.w)+df(h.w-1))/2*stretch
				delete(n.iset, h.w-1)
				heappush(&n.heap, heaptriple{mag(a.X[h.w-1]), h.w - 1, 0})
			}
		}

	}

	for i := range out {
		{
			p := G[`phasogram`].([][]float64)
			o := G[`origphase`].([][]float64)
			p[len(p)-1][i] = princarg(a.Frame[i])
			// o[len(p)-1][i] = o[len(p)-1][max(0, i-1)] + angle(a.P[i])
			o[len(p)-1][i] = angle(a.P[i])
		}
		out[i] = poltocar(mag(a.X[i]), a.Frame[i])
	}
	n.prevstretch = stretch
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
