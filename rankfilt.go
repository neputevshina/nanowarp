// Based on SciPy: scipy/ndimage/src/_rank_filter_1d.cpp
// Copyright (c) 2011 ashelly.myopenid.com under
// <http://www.opensource.org/licenses/mit-license>
// Modified in 2024 by Gideon Genadi Kogan
// Modified in 2025 by neputevshina
// Rewritten by ChatGPT to Go in 2025 for use in Nanowarp

package nanowarp

import (
	"golang.org/x/exp/constraints"
)

type mediatorMode int

const (
	NEAREST  mediatorMode = 0
	WRAP     mediatorMode = 1
	REFLECT  mediatorMode = 2
	MIRROR   mediatorMode = 3
	CONSTANT mediatorMode = 4
)

// Mediator - rank keeping structure
type Mediator[T constraints.Ordered] struct {
	pos   []int // index into `heap` for each value (signed position)
	heap  []int // underlying storage
	maxN  int   // allocated size
	N     int   // effective size (window length)
	idx   int   // position in circular queue
	minCt int   // count of items in min heap
	maxCt int   // count of items in max heap

	rawrank float64
	rank    int
	data    []T
}

// helper to access heap entries with signed index i (can be negative)
// maps i to heaporigin[i + rank]
func (m *Mediator[T]) getHeap(i int) int {
	return m.heap[i+m.rank]
}

func (m *Mediator[T]) setHeap(i int, val int) {
	m.heap[i+m.rank] = val
}

// returns true if data[heap[i]] < data[heap[j]]
func mmless[T constraints.Ordered](data []T, m *Mediator[T], i, j int) bool {
	return data[m.getHeap(i)] < data[m.getHeap(j)]
}

// swaps items i & j in heap, maintains pos
func mmexchange[T constraints.Ordered](m *Mediator[T], i, j int) {
	m.setHeap(i, m.getHeap(i))
	m.setHeap(j, m.getHeap(j))
	m.pos[m.getHeap(i)] = i
	m.pos[m.getHeap(j)] = j
}

// swaps items i & j if heap[i] < heap[j]; returns true if swapped
func mmCmpExch[T constraints.Ordered](data []T, m *Mediator[T], i, j int) bool {
	if mmless(data, m, i, j) {
		mmexchange(m, i, j)
		return true
	}
	return false
}

// maintains minheap property for all items below i.
func minSortDown[T constraints.Ordered](data []T, m *Mediator[T], i int) {
	for i *= 2; i <= m.minCt; i *= 2 {
		if i < m.minCt && mmless(data, m, i+1, i) {
			i++
		}
		if !mmCmpExch(data, m, i, i/2) {
			break
		}
	}
}

// maintains maxheap property for all items below i. (negative indexes)
func maxSortDown[T constraints.Ordered](data []T, m *Mediator[T], i int) {
	for i *= 2; i >= -m.maxCt; i *= 2 {
		if i > -m.maxCt && mmless(data, m, i, i-1) {
			i--
		}
		if !mmCmpExch(data, m, i/2, i) {
			break
		}
	}
}

// maintains minheap property for all items above i, including the rank
// returns true if rank changed (i reached 0)
func minSortUp[T constraints.Ordered](data []T, m *Mediator[T], i int) bool {
	for i > 0 && mmCmpExch(data, m, i, i/2) {
		i /= 2
	}
	return i == 0
}

// maintains maxheap property for all items above i, including rank
// returns true if rank changed (i reached 0)
func maxSortUp[T constraints.Ordered](data []T, m *Mediator[T], i int) bool {
	for i < 0 && mmCmpExch(data, m, i/2, i) {
		i /= 2
	}
	return i == 0
}

// Creates new Mediator: to calculate `nItems` running rank.
func MediatorNew[T constraints.Ordered](maxnItems, nItems int, rank float64) *Mediator[T] {
	if nItems > maxnItems {
		panic("nItems > maxnItems")
	}
	m := &Mediator[T]{}
	m.pos = make([]int, maxnItems)
	m.heap = make([]int, maxnItems)
	m.data = make([]T, maxnItems)
	m.maxN = maxnItems
	m.N = nItems
	m.rawrank = rank
	MediatorReset(m, nItems)
	return m
}

func MediatorReset[T constraints.Ordered](m *Mediator[T], nItems int) {
	m.N = nItems
	m.rank = int(m.rawrank * float64(nItems))
	if m.rank < 0 {
		m.rank = 0
	}
	if m.rank > m.maxN {
		m.rank = m.maxN
	}
	// heap virtual base is at index 0, mapped to heaporigin[rank]
	// initialize counts
	m.idx = 0
	m.minCt = m.N - m.rank - 1
	if m.minCt < 0 {
		m.minCt = 0
	}
	m.maxCt = m.rank
	// initialize pos and heap mapping
	for i := 0; i < nItems; i++ {
		p := i - m.rank
		m.pos[i] = p
		m.setHeap(p, i)
	}
	// zero data for first nItems
	var zero T
	for i := 0; i < nItems; i++ {
		m.data[i] = zero
	}
}

// MediatorDelete - not strictly necessary in Go (garbage collected) but provided to zero slices
func MediatorDelete[T constraints.Ordered](m *Mediator[T]) {
	if m == nil {
		return
	}
	m.heap = nil
	m.pos = nil
	m.data = nil
}

// Inserts item, maintains rank in O(lg nItems)
func MediatorInsert[T constraints.Ordered](data []T, m *Mediator[T], v T) {
	p := m.pos[m.idx]
	old := data[m.idx]
	data[m.idx] = v
	m.idx++
	if m.idx == m.N {
		m.idx = 0
	}

	if p > 0 { // new item is in minHeap
		if v > old {
			minSortDown(data, m, p)
			return
		}
		if minSortUp(data, m, p) && mmCmpExch(data, m, 0, -1) {
			maxSortDown(data, m, -1)
		}
	} else if p < 0 { // new item is in maxheap
		if v < old {
			maxSortDown(data, m, p)
			return
		}
		if maxSortUp(data, m, p) && mmCmpExch(data, m, 1, 0) {
			minSortDown(data, m, 1)
		}
	} else { // p == 0, new item is at rank
		if maxSortUp(data, m, -1) {
			maxSortDown(data, m, -1)
		}
		if minSortUp(data, m, 1) {
			minSortDown(data, m, 1)
		}
	}
}

// _rank_filter: applies a running rank filter over input slice `in` of length len,
// window size win, output to `out`. `mode` controls boundary handling, `cval` is constant value for CONSTANT mode, `origin` shifts the window center.
func _rank_filter[T constraints.Ordered](m *Mediator[T], in []T, length int, win int, out []T, mode mediatorMode, cval T, origin int) {
	if m == nil {
		panic("null Mediator")
	}
	if win <= 0 {
		panic("window must be > 0")
	}
	MediatorReset(m, win)
	data := m.data

	lim := (win-1)/2 - origin
	lim2 := length - lim
	var offset int

	switch mode {
	case REFLECT:
		for i := win - lim - 1; i > -1; i-- {
			MediatorInsert(data, m, in[i])
		}
	case CONSTANT:
		for i := win - lim; i > 0; i-- {
			MediatorInsert(data, m, cval)
		}
	case NEAREST:
		for i := win - lim; i > 0; i-- {
			MediatorInsert(data, m, in[0])
		}
	case MIRROR:
		for i := win - lim; i > 0; i-- {
			MediatorInsert(data, m, in[i])
		}
	case WRAP:
		if win%2 == 0 {
			offset = 2
		} else {
			offset = 0
		}
		start := length - lim - offset - 2*origin
		if start < 0 {
			start = 0
		}
		for i := start; i < length; i++ {
			MediatorInsert(data, m, in[i])
		}
	}

	for i := 0; i < lim; i++ {
		MediatorInsert(data, m, in[i])
	}
	for i := lim; i < length; i++ {
		MediatorInsert(data, m, in[i])
		out[i-lim] = data[m.getHeap(0)]
	}

	switch mode {
	case REFLECT:
		arrLenThresh := length - 1
		for i := 0; i < lim; i++ {
			MediatorInsert(data, m, in[arrLenThresh-i])
			out[lim2+i] = data[m.getHeap(0)]
		}
	case CONSTANT:
		for i := 0; i < lim; i++ {
			MediatorInsert(data, m, cval)
			out[lim2+i] = data[m.getHeap(0)]
		}
	case NEAREST:
		arrLenThresh := length - 1
		for i := 0; i < lim; i++ {
			MediatorInsert(data, m, in[arrLenThresh])
			out[lim2+i] = data[m.getHeap(0)]
		}
	case MIRROR:
		arrLenThresh := length - 2
		for i := 0; i < lim+1; i++ {
			MediatorInsert(data, m, in[arrLenThresh-i])
			out[lim2+i] = data[m.getHeap(0)]
		}
	case WRAP:
		for i := 0; i < lim; i++ {
			MediatorInsert(data, m, in[i])
			out[lim2+i] = data[m.getHeap(0)]
		}
	}
}

// // Example usage (for reference)
// func main() {
// 	// Example with int type
// 	in := []int{5, 2, 8, 3, 7, 1, 4}
// 	length := len(in)
// 	win := 3
// 	out := make([]int, length)
// 	m := MediatorNew[int](20, win, 0.5) // rank 0.5 => median-ish
// 	_rank_filter[int](m, in, length, win, out, NEAREST, 0, 0)
// 	fmt.Println("input:", in)
// 	fmt.Println("output:", out)
// }
