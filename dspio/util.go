package dspio

import (
	"fmt"
	"os"

	"golang.org/x/exp/constraints"
)

func println(v ...any) { fmt.Fprintln(os.Stderr, v...) }

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

func copy2(dst [][]float64, src [][]float64) int {
	if len(dst) != len(src) {
		panic(`copy2: channel count mismatch`)
	}
	n := copy(dst[0], src[0])
	for ch := range src[1:] {
		nn := copy(dst[ch+1], src[ch+1])
		if n != nn {
			panic(`copy2: buffer size mismatch`)
		}
	}
	return n
}

// func dump(name string, data []float64, fs int) {
// 	file, err := os.Create(name)
// 	defer file.Close()
// 	if err != nil {
// 		panic(err)
// 	}

// 	wr := wav.NewWriter(file, uint32(len(data)), 1, uint32(fs), 32, true)
// 	nbuf := 2048
// 	buf := make([]wav.Sample, 0, nbuf)
// 	for i := 0; i < len(data); i += nbuf {
// 		buf = buf[:0]
// 		for j := i; j < min(i+nbuf, len(data)); j++ {
// 			buf = append(buf, wav.Sample{Values: [2]int{
// 				int(math.Float32bits(float32(data[j])))}})
// 		}
// 		err := wr.WriteSamples(buf)
// 		if err != nil {
// 			panic(err)
// 		}
// 	}
// }

func boolfloat(b bool) float64 {
	if b {
		return 1
	}
	return 0
}
