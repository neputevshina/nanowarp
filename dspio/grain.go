package dspio

import "fmt"

type GrainReader struct {
	Hop       int // Can be freely changed outside
	n         int
	nch       int
	inbuss    [][]float64
	slicebuss [][]float64
	ar        SignalReader
}

func NewGrainReader(nfft, hop int, ar SignalReader) (g *GrainReader) {
	nch := ar.NchRead()
	g = &GrainReader{
		Hop:       hop,
		n:         nfft,
		nch:       nch,
		inbuss:    make([][]float64, nch),
		slicebuss: make([][]float64, nch),
		ar:        ar}
	for ch := range nch {
		g.inbuss[ch] = make([]float64, nfft)
	}
	return
}

// SignalRead reads next overlapped grains from input.
//
// SignalRead expects grains[ch] for every ch to contain at least nfft samples.
func (g *GrainReader) SignalRead(prr error, grain [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(grain) != g.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, g.nch, len(grain)))
	}
	for ch := range g.inbuss {
		copy(g.inbuss[ch], g.inbuss[ch][g.Hop:])
		g.slicebuss[ch] = g.inbuss[ch][g.n-g.Hop:]
	}
	s := 0
	for s < g.Hop {
		n, err = g.ar.SignalRead(nil, g.slicebuss[s:])
		if err != nil {
			return s, err
		}
		s += n
	}
	for ch := range g.inbuss {
		copy(grain[ch], g.inbuss[ch])
	}
	return s, nil
}

func (g *GrainReader) NchRead() int { return g.nch }

var _ SignalReader = &GrainReader{}

type GrainWriter struct {
	Hop       int // Can be freely changed outside
	n         int
	nch       int
	outbuss   [][]float64
	slicebuss [][]float64
	aw        SignalWriter
}

func NewGrainWriter(nfft, hop int, aw SignalWriter) (g *GrainWriter) {
	nch := aw.NchWrite()
	g = &GrainWriter{
		Hop:       hop,
		n:         nfft,
		nch:       nch,
		outbuss:   make([][]float64, nch),
		slicebuss: make([][]float64, nch),
		aw:        aw}
	for ch := range nch {
		g.outbuss[ch] = make([]float64, nfft*2)
	}
	return
}

// SignalWrite writes next overlap-added grains to the output.
//
// SignalWrite takes no more than nfft samples from each channel of grains.
func (g *GrainWriter) SignalWrite(prr error, grain [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(grain) != g.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, g.nch, len(grain)))
	}
	for ch := range g.outbuss {
		add(g.outbuss[ch][g.Hop:], grain[ch])
		g.slicebuss[ch] = g.outbuss[ch][:g.Hop]
	}

	s := 0
	for s < g.Hop {
		n, err = g.aw.SignalWrite(nil, g.slicebuss[s:])
		if err != nil {
			return s, err
		}
		s += n
	}

	for ch := range g.outbuss {
		copy(g.outbuss[ch], g.outbuss[ch][g.Hop:])
	}
	return s, nil
}

func (g *GrainWriter) NchWrite() int { return g.nch }

var _ SignalWriter = &GrainWriter{}
