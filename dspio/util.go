package dspio

import (
	"fmt"
	"math"
	"os"

	"github.com/neputevshina/nanowarp/wav"
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

func dump(name string, data []float64, fs int) {
	file, err := os.Create(name)
	defer file.Close()
	if err != nil {
		panic(err)
	}

	wr := wav.NewWriter(file, uint32(len(data)), 1, uint32(fs), 32, true)
	nbuf := 2048
	buf := make([]wav.Sample, 0, nbuf)
	for i := 0; i < len(data); i += nbuf {
		buf = buf[:0]
		for j := i; j < min(i+nbuf, len(data)); j++ {
			buf = append(buf, wav.Sample{Values: [2]int{
				int(math.Float32bits(float32(data[j])))}})
		}
		err := wr.WriteSamples(buf)
		if err != nil {
			panic(err)
		}
	}
}
