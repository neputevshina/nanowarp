package nanowarp

import (
	"container/heap"
	"fmt"
	"math"
	"math/cmplx"
	"os"
)

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

func (n *Nanowarp) Process1(in []float64, out []float64, stretch float64) {
	a := &n.a
	ih := int(math.Floor(float64(n.hop) / stretch))
	oh := int(math.Floor(float64(n.hop)))
	fmt.Fprintln(os.Stderr, `inhop:`, ih, `nbuf:`, n.nbuf)
	fmt.Fprintln(os.Stderr, `outhop:`, oh)

	adv := 0
	for i := 0; i < len(in); i += ih {
		enfft := func(i int, x []complex128, w []float64) {
			zero(a.S)
			copy(a.S, in[i:i+n.nbuf])
			mul(a.S, w)
			n.fft.Coefficients(x, a.S)
		}

		enfft(i, a.X, a.W)
		enfft(i, a.Xd, a.Wd)
		enfft(i, a.Xt, a.Wt)

		if i == 0 {
			copy(a.P, a.X)
			continue
		}

		// Begin of PGHI

		a := &n.a
		n.heap = make(hp, n.nbins)
		clear(n.arm)

		for j := range a.X {
			n.arm[j] = true
			n.heap[j] = heaptriple{mag(a.P[j]), j, -1}
		}
		heap.Init(&n.heap)

		olap := float64(n.nbuf / n.hop)
		padv := func(j int) float64 {
			// This phase advance formula is a product of 10 hours of trial and error.
			return (math.Pi*float64(j) + imag(a.Xd[j]/a.X[j])) / olap
		}
		fadv := func(j int) float64 {
			return (-math.Pi - real(a.Xt[j]/a.X[j])/float64(len(a.X))*2*math.Pi) / 2
		}

		for len(n.heap) > 0 {

			h := heap.Pop(&n.heap).(heaptriple)
			w := h.w
			switch h.t {
			case -1:
				if n.arm[w] {
					a.Phase[w] = princarg(cmplx.Phase(a.P[w]) + padv(w))
					n.arm[w] = false
					heap.Push(&n.heap, heaptriple{mag(izero(a.X, w)), w, 0})
				}
			case 0:
				if w > 1 && n.arm[w-1] {
					a.Phase[w-1] = princarg(a.Phase[w] - fadv(w))
					n.arm[w-1] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w-1]), w - 1, 0})
				}
				if w < n.nbins-1 && n.arm[w+1] {
					a.Phase[w+1] = princarg(a.Phase[w] + fadv(w))
					n.arm[w+1] = false
					heap.Push(&n.heap, heaptriple{mag(a.X[w+1]), w + 1, 0})
				}
			}
		}

		for j := range a.Phase {
			a.A[j] = cmplx.Rect(mag(a.X[j]), a.Phase[j])
		}
		copy(a.P, a.A)

		// End of PGHI

		n.fft.Sequence(a.S, a.P)
		for j := range a.S {
			a.S[j] /= n.norm
			a.S[j] /= stretch
		}
		mul(a.S, a.W)
		add(out[adv:adv+n.nbuf], a.S)
		adv += oh
	}
}
