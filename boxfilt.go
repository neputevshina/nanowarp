package nanowarp

type boxfilt struct {
	ring []float64
	n, j int
	avg  float64
}

func boxfiltNew(nItems int) *boxfilt {
	f := &boxfilt{}
	f.n = nItems
	f.ring = make([]float64, 2*nItems)
	f.j = 0
	return f
}

func (f *boxfilt) Insert(v float64) {
	conj := func(j int) int { return (j + f.n) % f.n }
	f.avg += (v - f.ring[f.j]) / float64(f.n)
	f.ring[f.j] = v
	f.ring[conj(f.j)] = v
	f.j++
	f.j = conj(f.j)
}

func (f *boxfilt) Take() (avg, delay float64) {
	return f.avg, f.ring[f.j]
}
