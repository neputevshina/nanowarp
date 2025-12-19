package nanowarp

// This is a manually specialized copy of container/heap.

type heaptriple struct {
	mag  float64
	w, t int
}
type hp []heaptriple

func (h hp) Len() int           { return len(h) }
func (h hp) Less(i, j int) bool { return h[i].mag > h[j].mag }
func (h hp) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *hp) Push(x any) {
	*h = append(*h, x.(heaptriple))
}
func (h *hp) Pop() any {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

// heapInit establishes the heap invariants required by the other routines in this package.
// heapInit is idempotent with respect to the heap invariants
// and may be called whenever the heap invariants may have been invalidated.
// The complexity is O(n) where n = h.Len().
func heapInit(h *hp) {
	// heapify
	n := h.Len()
	for i := n/2 - 1; i >= 0; i-- {
		heapdown(h, i, n)
	}
}

// heapPush pushes the element x onto the heap.
// The complexity is O(log n) where n = h.Len().
func heapPush(h *hp, x any) {
	h.Push(x)
	heapup(h, h.Len()-1)
}

// heapPop removes and returns the minimum element (according to Less) from the heap.
// The complexity is O(log n) where n = h.Len().
// heapPop is equivalent to [heapRemove](h, 0).
func heapPop(h *hp) any {
	n := h.Len() - 1
	h.Swap(0, n)
	heapdown(h, 0, n)
	return h.Pop()
}

// heapRemove removes and returns the element at index i from the heap.
// The complexity is O(log n) where n = h.Len().
func heapRemove(h *hp, i int) any {
	n := h.Len() - 1
	if n != i {
		h.Swap(i, n)
		if !heapdown(h, i, n) {
			heapup(h, i)
		}
	}
	return h.Pop()
}

// heapFix re-establishes the heap ordering after the element at index i has changed its value.
// Changing the value of the element at index i and then calling heapFix is equivalent to,
// but less expensive than, calling [heapRemove](h, i) followed by a Push of the new value.
// The complexity is O(log n) where n = h.Len().
func heapFix(h *hp, i int) {
	if !heapdown(h, i, h.Len()) {
		heapup(h, i)
	}
}

func heapup(h *hp, j int) {
	for {
		i := (j - 1) / 2 // parent
		if i == j || !h.Less(j, i) {
			break
		}
		h.Swap(i, j)
		j = i
	}
}

func heapdown(h *hp, i0, n int) bool {
	i := i0
	for {
		j1 := 2*i + 1
		if j1 >= n || j1 < 0 { // j1 < 0 after int overflow
			break
		}
		j := j1 // left child
		if j2 := j1 + 1; j2 < n && h.Less(j2, j1) {
			j = j2 // = 2*i + 2  // right child
		}
		if !h.Less(j, i) {
			break
		}
		h.Swap(i, j)
		i = j
	}
	return i > i0
}
