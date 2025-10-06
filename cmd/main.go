package main

import (
	"io"
	"math"
	"os"

	"github.com/neputevshina/nanowarp"
	"github.com/youpy/go-wav"
)

func main() {
	file, _ := os.Open(`ticktock.wav`)
	rd := wav.NewReader(file)
	data := []float64{}
	for {
		samples, err := rd.ReadSamples()
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			data = append(data, rd.FloatValue(sample, 0))
		}
	}
	nw := nanowarp.New()
	n := 2.
	out := make([]float64, int(float64(len(data))*n)+12000)
	nw.Process(data, out, n)

	for i := range out {
		out[i] = out[i] / math.Abs(out[i]) * max(0, min(math.Abs(out[i]), 1))
	}

	file, err := os.Create(`ticktock-x2.wav`)
	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(out)), 1, 48000, 16)
	for _, sa := range out {
		err := wr.WriteSamples([]wav.Sample{{Values: [2]int{int(sa * math.Pow(2, 16-1))}}})
		if err != nil {
			panic(err)
		}
	}
	file.Close()

}
