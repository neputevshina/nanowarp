package oscope

import (
	"container/heap"
	"fmt"
	"math"
	"path"
	"runtime"
	"slices"
	"sort"
)

var dists = map[where]*dist{}

type histopair struct {
	v float64
	n int
}
type histoheap []histopair

func (h *histoheap) Len() int               { return len(*h) }
func (h *histoheap) Less(i int, j int) bool { return (*h)[i].v > (*h)[j].v }
func (h *histoheap) Swap(i int, j int) {
	a, b := (*h)[i], (*h)[j]
	(*h)[i], (*h)[j] = b, a
}
func (h *histoheap) Push(x any) { *h = append(*h, x.(histopair)) }
func (h *histoheap) Pop() (v any) {
	v = (*h)[len(*h)-1]
	*h = (*h)[:len(*h)-1]
	return v
}

var _ heap.Interface = &histoheap{}

type dist struct {
	outfile string
	heap    histoheap
}

func Histo(v float64, etc ...Etc) {
	if !Enable {
		return
	}
	pc, file, line, ok := runtime.Caller(1)
	wh := [2]uintptr{pc, getGID()}
	d, k := dists[wh]
	if !k {
		if !ok {
			panic(`oscope.Histo: can't identify a function that is not traceable on the stack`)
		}
		fn := ""
		_, e, ok := findetc[Name](etc)
		if ok {
			fn = string(e)
		} else {
			fn = fmt.Sprintf("%s:%d(%d)", path.Base(file), line, wh[1])
		}

		dists[wh] = &dist{
			outfile: fn,
			heap:    make([]histopair, 0),
		}
		d = dists[wh]
	}

	// Quantize v to set bits of mantissa, 10 by default.
	_, q, ok := findetc[Quantization](etc)
	if !ok {
		q = 10
	}
	vq := math.Float64frombits(math.Float64bits(v) >> (53 - q) << (53 - q))

	i := slices.IndexFunc(d.heap, func(p histopair) bool { return p.v == vq })
	if i >= 0 {
		d.heap[i].n++
	} else {
		heap.Push(&d.heap, histopair{v: vq, n: 1})
	}
}

func dumpDist(_ error, w *dist, _ string) error {
	fmt.Println(`Distribution for `, w.outfile)
	sort.Sort(&w.heap)
	for _, p := range w.heap {
		fmt.Println(p.v, "\t", p.n)
	}
	return nil
}
