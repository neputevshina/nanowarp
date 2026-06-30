package nanowarp

import (
	"gonum.org/v1/gonum/dsp/fourier"
)

// Masuyama, Yoshiki, Kohei Yatabe, and Yasuhiro Oikawa.
// "Griffin–Lim like phase recovery via alternating direction method of multipliers."
// IEEE Signal Processing Letters 26.1 (2018): 184-188.
// https://ieeexplore.ieee.org/iel7/97/4358004/08552369.pdf

type admm struct {
	fft     *fourier.FFT
	nbuf    int
	s, sbig []float64
	norm    float64
	hop     int
	Wf, Wr  []float64

	admmbufs
}

type admmbufs struct { // DELETEME
	Z, X, U, Y, Yp [][]complex128 `size:"lah"`
}

type stftTensor struct { // TODO(neputevshina): USE ME
	data    []complex128
	shape   [2]int
	padding int
	view    [][]complex128
}

func (a *admm) stft(dst [][]complex128, sa []float64) [][]complex128 {
	nbuf := a.nbuf
	for j := range dst {
		i := j * a.hop
		copy(a.s, sa[i:i+nbuf])
		mul(a.s, a.Wf)
		a.fft.Coefficients(dst[j], a.s)
	}
	return dst
}

func (a *admm) istft(dst []float64, fr [][]complex128) []float64 {
	nbuf := a.nbuf
	for j := range fr {
		i := j * a.hop
		a.fft.Sequence(a.s, fr[j])
		mul(a.s, a.Wr)
		scale(a.s, 1/a.norm)
		add(dst[i:i+nbuf], a.s)
	}
	return dst
}

func (a *admm) project(dst [][]complex128) [][]complex128 {
	clear(a.sbig)
	return a.stft(dst, a.istft(a.sbig[:a.hop*8+2*a.nbuf], dst))
}

// gla implements basic Griffin-Lim phase reconstruction algorithm.
func (a *admm) gla(X [][]complex128, M [][]float64, known []bool, iterations int) {
	if a.U == nil { // TODO(neputevshina): create newGla and move to it in production nanowarp
		a.Z = make2[complex128](len(X), len(X[0]))
	}

	Z := a.Z
	for n := range X {
		copy(Z[n], X[n])
	}
	for range iterations {
		for i := range Z {
			for j := range Z[0] {
				X[i][j] = complex(M[i][j], 0) * norm(X[i][j])
			}
		}
		X = a.project(X)

		for i := range Z {
			for j := range Z[0] {
				if known[i] {
					X[i][j] = Z[i][j]
				}
			}
		}
	}
	// waveform.Dump(nil, a.sbig[:a.hop*8+2*a.nbuf])
}

// Based on ADMMGLA.m from original paper's repo.
func (a *admm) admm(X [][]complex128, M [][]float64, known []bool, iterations int, ρ float64) {
	if a.U == nil { // TODO(neputevshina): create newAdmm and move to it in production nanowarp
		a.Z = make2[complex128](len(X), len(X[0]))
		a.U = make2[complex128](len(X), len(X[0]))
		a.Y = make2[complex128](len(X), len(X[0]))
		a.Yp = make2[complex128](len(X), len(X[0]))
	}

	Z, U, Y, Yp := a.Z, a.U, a.Y, a.Yp
	for n := range X {
		copy(Z[n], X[n])
		// copy(X[n],X[n])
		clear(U[n])
	}
	for range iterations {
		for i := range Z {
			for j := range Z[0] {
				if !known[i] { // Retain the phase if known, |X| = M
					X[i][j] = complex(M[i][j], 0) * norm(Z[i][j]-U[i][j]) // Pc2
				} else {
					Z[i][j] = X[i][j]
				}
				Y[i][j] = X[i][j] + U[i][j]
				Yp[i][j] = Y[i][j]
			}
		}
		Yp = a.project(Yp)
		// waveform.Dump(nil, a.sbig[:a.hop*8+2*a.nbuf])

		for i := range Z {
			for j := range Z[0] {
				ρ := complex(ρ, 0)
				Z[i][j] = (ρ*Y[i][j] + Yp[i][j]) / (1 + ρ)
				U[i][j] += X[i][j] - Z[i][j]
			}
		}
	}
}
