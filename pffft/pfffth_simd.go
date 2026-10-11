//go:build goexperiment.simd

package pffft

type pffftSetup struct {
	n int
	// nb of complex simd vectors (N/4 if PFFFT_COMPLEX, N/8 if PFFFT_REAL)
	// hold the decomposition into small integers of N
	ncvec int
	// N , number of factors, factors (admitted values: 2, 3, 4 ou 5)
	ifac [IfacMaxSize]int
	transform
	data    []lane    // allocated room for twiddle coefs
	e       []float32 // points into 'data' , N/4*3 elements
	twiddle []float32 // points into 'data', N/4 elements
}

type transform uint

const (
	transformReal transform = iota
	transformComplex
)
