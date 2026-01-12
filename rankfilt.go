// Based on SciPy: scipy/ndimage/src/_rank_filter_1d.cpp
// Copyright (c) 2011 ashelly.myopenid.com under
// <http://www.opensource.org/licenses/mit-license>
// Modified in 2024 by Gideon Genadi Kogan
// Modified in 2025 by neputevshina
// Rewritten by neputevshina and ChatGPT to Go in 2025 for use in Nanowarp

package nanowarp

import "golang.org/x/exp/constraints"

type mediatorMode int

const (
	mNEAREST  mediatorMode = 0
	mWRAP     mediatorMode = 1
	mREFLECT  mediatorMode = 2
	mMIRROR   mediatorMode = 3
	mCONSTANT mediatorMode = 4
)

// mediator - rank keeping structure
type mediator[T constraints.Ordered, I any] struct {
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
	more    []I
}

// Creates new Mediator: to calculate `nItems` running rank.
func mediatorNew[T constraints.Ordered, I any](maxnItems, nItems int, rank float64) *mediator[T, I] {
	if nItems > maxnItems {
		panic("nItems > maxnItems")
	}
	m := &mediator[T, I]{}
	m.pos = make([]int, maxnItems*2)
	m.heap = make([]int, maxnItems*2)
	m.data = make([]T, maxnItems*2)
	m.more = make([]I, maxnItems*2)
	m.maxN = maxnItems
	m.N = nItems
	m.rawrank = rank
	m.Reset(nItems)
	return m
}

// helper to access heap entries with signed index i (can be negative)
// maps i to heaporigin[i + rank]
func (m *mediator[T, I]) getHeap(i int) int {
	return m.heap[i+m.rank]
}

func (m *mediator[T, I]) setHeap(i int, val int) {
	m.heap[i+m.rank] = val
}

// returns true if data[heap[i]] < data[heap[j]]
func (m *mediator[T, I]) less(i, j int) bool {
	// return data[m.getHeap(i)] < data[m.getHeap(j)]
	return m.data[m.getHeap(i)] < m.data[m.getHeap(j)]
}

// swaps items i & j in heap, maintains pos
func (m *mediator[T, I]) xch(i, j int) {
	a, b := m.getHeap(i), m.getHeap(j)
	m.setHeap(i, b)
	m.setHeap(j, a)
	m.pos[m.getHeap(i)] = i
	m.pos[m.getHeap(j)] = j
}

// swaps items i & j if heap[i] < heap[j]; returns true if swapped
func (m *mediator[T, I]) lessxch(i, j int) bool {
	if m.less(i, j) {
		m.xch(i, j)
		return true
	}
	return false
}

// maintains minheap property for all items below i.
func (m *mediator[T, I]) minSortDown(i int) {
	for i *= 2; i <= m.minCt; i *= 2 {
		if i < m.minCt && m.less(i+1, i) {
			i++
		}
		if !m.lessxch(i, i/2) {
			break
		}
	}
}

// maintains maxheap property for all items below i. (negative indexes)
func (m *mediator[T, I]) maxSortDown(i int) {
	for i *= 2; i >= -m.maxCt; i *= 2 {
		if i > -m.maxCt && m.less(i, i-1) {
			i--
		}
		if !m.lessxch(i/2, i) {
			break
		}
	}
}

// maintains minheap property for all items above i, including the rank
// returns true if rank changed (i reached 0)
func (m *mediator[T, I]) minSortUp(i int) bool {
	for i > 0 && m.lessxch(i, i/2) {
		i /= 2
	}
	return i == 0
}

// maintains maxheap property for all items above i, including rank
// returns true if rank changed (i reached 0)
func (m *mediator[T, I]) maxSortUp(i int) bool {
	for i < 0 && m.lessxch(i/2, i) {
		i /= 2
	}
	return i == 0
}

func (m *mediator[T, I]) Reset(nItems int) {
	m.N = nItems
	m.rank = int(m.rawrank * float64(nItems))
	// heap virtual base is at index 0, mapped to heap[rank]
	// initialize counts
	m.idx = 0
	m.minCt = m.N - m.rank - 1
	m.maxCt = m.rank
	for i := range nItems {
		p := i - m.rank
		m.pos[i] = p
		m.setHeap(p, i)
	}
	var zero T
	for i := range nItems {
		m.data[i] = zero
	}
}

func (m *mediator[T, I]) Take() (v T, a I) {
	c := m.heap[m.rank]
	return m.data[c], m.more[c]
}

// Inserts item, maintains rank in O(lg nItems)
func (m *mediator[T, I]) Insert(v T, a I) {
	p := m.pos[m.idx]
	old := m.data[m.idx]
	m.data[m.idx] = v
	m.more[m.idx] = a
	m.idx++
	if m.idx == m.N {
		m.idx = 0
	}

	if p > 0 { // new item is in minHeap
		if v > old {
			m.minSortDown(p)
			return
		}
		if m.minSortUp(p) && m.lessxch(0, -1) {
			m.maxSortDown(-1)
		}
	} else if p < 0 { // new item is in maxheap
		if v < old {
			m.maxSortDown(p)
			return
		}
		if m.maxSortUp(p) && m.lessxch(1, 0) {
			m.minSortDown(1)
		}
	} else { // p == 0, new item is at rank
		if m.maxSortUp(-1) {
			m.maxSortDown(-1)
		}
		if m.minSortUp(1) {
			m.minSortDown(1)
		}
	}
}

// filt applies a running rank filter over input slice `in` of length len,
// window size win, output to `out`. `mode` controls boundary handling, `cval` is a
// constant value for CONSTANT mode, `origin` shifts the window center.
func (m *mediator[T, I]) filt(in []T, win int, out []T, mode mediatorMode, cval T, origin int) {
	if m == nil {
		panic("null Mediator")
	}
	if win <= 0 {
		panic("window must be > 0")
	}
	m.Reset(win)
	data := m.data
	length := len(in)

	lim := (win-1)/2 - origin
	lim2 := length - lim
	var offset int

	var z I

	switch mode {
	case mREFLECT:
		for i := win - lim - 1; i > -1; i-- {
			m.Insert(in[i], z)
		}
	case mCONSTANT:
		for i := win - lim; i > 0; i-- {
			m.Insert(cval, z)
		}
	case mNEAREST:
		for i := win - lim; i > 0; i-- {
			m.Insert(in[0], z)
		}
	case mMIRROR:
		for i := win - lim; i > 0; i-- {
			m.Insert(in[i], z)
		}
	case mWRAP:
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
			m.Insert(in[i], z)
		}
	}

	for i := 0; i < lim; i++ {
		m.Insert(in[i], z)
	}
	for i := lim; i < length; i++ {
		m.Insert(in[i], z)
		out[i-lim] = data[m.getHeap(0)]
	}

	switch mode {
	case mREFLECT:
		arrLenThresh := length - 1
		for i := 0; i < lim; i++ {
			m.Insert(in[arrLenThresh-i], z)
			out[lim2+i] = data[m.getHeap(0)]
		}
	case mCONSTANT:
		for i := 0; i < lim; i++ {
			m.Insert(cval, z)
			out[lim2+i] = data[m.getHeap(0)]
		}
	case mNEAREST:
		arrLenThresh := length - 1
		for i := 0; i < lim; i++ {
			m.Insert(in[arrLenThresh], z)
			out[lim2+i] = data[m.getHeap(0)]
		}
	case mMIRROR:
		arrLenThresh := length - 2
		for i := 0; i < lim+1; i++ {
			m.Insert(in[arrLenThresh-i], z)
			out[lim2+i] = data[m.getHeap(0)]
		}
	case mWRAP:
		for i := 0; i < lim; i++ {
			m.Insert(in[i], z)
			out[lim2+i] = data[m.getHeap(0)]
		}
	}
}
