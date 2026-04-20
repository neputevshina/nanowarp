package nanowarp

type AudioReader interface {
	NchRead() int
	AudioRead(prr error, buf [][]float64) (n int, err error)
}

type AudioWriter interface {
	NchWrite() int
	AudioWrite(prr error, buf [][]float64) (n int, err error)
}

type GrainReader struct {
	Hop       int // Can be freely changed outside
	n         int
	inbuss    [][]float64
	slicebuss [][]float64
	ar        AudioReader
}

func NewGrainReader(nfft, hop int, ar AudioReader) (g *GrainReader) {
	nch := ar.NchRead()
	g = &GrainReader{
		Hop:       hop,
		n:         nfft,
		inbuss:    make([][]float64, nch),
		slicebuss: make([][]float64, nch),
		ar:        ar}
	for ch := range nch {
		g.inbuss[ch] = make([]float64, nfft)
	}
	return
}

func (g *GrainReader) AudioRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	for ch := range g.inbuss {
		copy(g.inbuss[ch], g.inbuss[ch][g.Hop:])
		g.slicebuss[ch] = g.inbuss[ch][g.n-g.Hop:]
	}
	s := 0
	for s < g.Hop {
		n, err = g.ar.AudioRead(nil, g.slicebuss[s:])
		if err != nil {
			return s, err
		}
		s += n
	}
	return s, nil
}

func (g *GrainReader) NchRead() int {
	return g.ar.NchRead()
}

var _ AudioReader = &GrainReader{}

type GrainWriter struct {
}
