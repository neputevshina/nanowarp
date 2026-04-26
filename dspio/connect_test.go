package dspio

import (
	"math/rand"
	"testing"
)

func TestPipeBusy(t *testing.T) {
	w, r := Pipe(1, 8192)

}

type noiseOsc struct {
	i, p int
}

func (p *noiseOsc) NchRead() int { return 1 }

func (p *noiseOsc) SignalRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	for i := range buf[0] {
		buf[0][i] = rand.Float64()
	}
	return len(buf[0]), nil
}

var _ SignalReader = &noiseOsc{}
