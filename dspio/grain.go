package dspio

import "fmt"

type GrainReader struct {
	Hop       int // Can be freely changed outside
	n         int
	nch       int
	inbuss    [][]float64
	slicebuss [][]float64
	r         SignalReader
}

func NewGrainReader(nfft, hop int, r SignalReader) (g *GrainReader) {
	nch := r.NchRead()
	g = &GrainReader{
		Hop:       hop,
		n:         nfft,
		nch:       nch,
		inbuss:    make([][]float64, nch),
		slicebuss: make([][]float64, nch),
		r:         r}
	for ch := range nch {
		g.inbuss[ch] = make([]float64, nfft)
	}
	return

}

func (r *GrainReader) N() int { return r.n }

// SignalRead reads the next overlapped grain from input.
// It returns the amount of samples actually read from the input (at least r.Hop),
// but fills the grain with r.N()-r.Hop previous samples.
//
// SignalRead expects grain[ch] for every ch to contain at least nfft samples.
func (r *GrainReader) SignalRead(prr error, grain [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(grain) != r.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, r.nch, len(grain)))
	}
	for ch := range r.inbuss {
		copy(r.inbuss[ch], r.inbuss[ch][r.Hop:])
	}
	s := 0
	for s < r.Hop {
		for ch := range r.slicebuss {
			r.slicebuss[ch] = r.inbuss[ch][r.n-r.Hop:][s:]
		}
		n, err = r.r.SignalRead(nil, r.slicebuss)
		if err != nil {
			return s, err
		}
		s += n
	}
	for ch := range r.inbuss {
		copy(grain[ch], r.inbuss[ch])
	}
	return s, nil
}

func (r *GrainReader) NchRead() int { return r.nch }

var _ SignalReader = &GrainReader{}

type GrainWriter struct {
	Hop       int // Can be freely changed outside
	n         int
	nch       int
	outbuss   [][]float64
	slicebuss [][]float64
	w         SignalWriter
}

func NewGrainWriter(nfft, hop int, w SignalWriter) (g *GrainWriter) {
	nch := w.NchWrite()
	g = &GrainWriter{
		Hop:       hop,
		n:         nfft,
		nch:       nch,
		outbuss:   make([][]float64, nch),
		slicebuss: make([][]float64, nch),
		w:         w}
	for ch := range nch {
		g.outbuss[ch] = make([]float64, nfft*2)
	}
	return
}

func (w *GrainWriter) N() int { return w.n }

// SignalWrite writes the next overlap-added grain to the output.
// It returns the amount of samples written to the output (at least r.Hop),
// but adds the whole buffer to the accumulator, so data in grain is expected to be
// different each time.
//
// SignalWrite takes no more than nfft samples from each channel of grain.
func (w *GrainWriter) SignalWrite(prr error, grain [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(grain) != w.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, w.nch, len(grain)))
	}
	for ch := range w.outbuss {
		add(w.outbuss[ch][w.Hop:], grain[ch])
	}

	s := 0
	for s < w.Hop {
		for ch := range w.slicebuss {
			w.slicebuss[ch] = w.outbuss[ch][:w.Hop][s:]
		}
		n, err = w.w.SignalWrite(nil, w.slicebuss)
		if err != nil {
			return s, err
		}
		s += n
	}

	for ch := range w.outbuss {
		// Note that len(w.outbuss[ch]) == nfft*2.
		// Even if w.Hop is equal to nfft, which is the maximum allowable,
		// we still have the rest of the buffer padded with zeros.
		//
		// For typical grain sizes and channel configurations (nfft ∈ [256:16384],
		// stereo) the 2x memory overhead is not significant.
		// In most cases it will be in range from 512 (cross-synthesis) to 4096 (phase
		// vocoder-based pitch change/time scaling).
		// The most probable reason for worst case (262kB) is either linear-phase EQ or
		// spectrum visualization, both of which are generally operations that are not
		// composed with other processes in a pipeline, so overhead won't scale.
		copy(w.outbuss[ch], w.outbuss[ch][w.Hop:])
	}
	return s, nil
}

func (w *GrainWriter) NchWrite() int { return w.nch }

var _ SignalWriter = &GrainWriter{}
