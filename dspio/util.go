package dspio

import (
	"golang.org/x/exp/constraints"
)

func add[T constraints.Float](dst, src []T) {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] += src[i]
	}
}

func make2(nch, n int) [][]float64 {
	g := make([][]float64, nch)
	for ch := range g {
		g[ch] = make([]float64, n)
	}
	return g
}
