package nanowarp

import "gonum.org/v1/gonum/dsp/fourier"

// Masuyama, Yoshiki, Kohei Yatabe, and Yasuhiro Oikawa.
// "Griffin–Lim like phase recovery via alternating direction method of multipliers."
// IEEE Signal Processing Letters 26.1 (2018): 184-188.
// https://ieeexplore.ieee.org/iel7/97/4358004/08552369.pdf

type admm struct {
	fft     *fourier.FFT
	s, sbig []float64
	norm    float64
	hop     int
	Wf, Wr  []float64
}

type stftTensor struct {
	data    []complex128
	shape   [2]int
	padding int
	view    [][]complex128
}

func (a *admm) stft(dst [][]complex128, sa []float64) [][]complex128 {
	nfft := a.fft.Len()
	for i := 0; i < len(sa); i += a.hop {
		copy(a.s[max(0, nfft/2-i):], sa[max(0, i-nfft/2):min(len(sa), i+nfft/2)])
		mul(a.s, a.Wf)
		a.fft.Coefficients(dst[i/a.hop], a.s)
	}
	return dst
}

func (a *admm) istft(dst []float64, fr [][]complex128) []float64 {
	nfft := a.fft.Len()
	for j := range fr {
		i := j * a.hop
		a.fft.Sequence(a.s, fr[j])
		mul(a.s, a.Wr)
		scale(a.s, 1/a.norm)
		add(dst[max(0, i-nfft/2):min(len(dst), i+nfft/2)], a.s)
	}
	return dst
}

func (a *admm) project(dst [][]complex128) [][]complex128 {
	return a.stft(dst, a.istft(a.sbig, dst))
}

func (a *admm) admm(srcdst [][]complex128, M [][]float64, iterations int, ρ float64) [][]complex128 {
	Z := make2[complex128](len(srcdst), len(srcdst[0]))
	for n := range srcdst {
		copy(Z[n], srcdst[n])
	}
	U := make2[complex128](len(srcdst), len(srcdst[0]))
	Y := make2[complex128](len(srcdst), len(srcdst[0]))
	W := make2[complex128](len(srcdst), len(srcdst[0]))
	X := srcdst
	for range iterations {
		for i := range Z {
			for j := range Z[0] {
				// Pc2
				X[i][j] = complex(M[i][j], 0) * norm(Z[i][j]-U[i][j])
				Y[i][j] = X[i][j] + U[i][j]
				W[i][j] = Y[i][j]
			}
		}
		W = a.project(W)
		for i := range Z {
			for j := range Z[0] {
				ρ := complex(ρ, 0)
				Z[i][j] = (ρ*Y[i][j] + W[i][j]) / (1 + ρ)
				U[i][j] = U[i][j] + X[i][j] - Z[i][j]
			}
		}
	}
	return srcdst
}
