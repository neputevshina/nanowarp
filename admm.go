package nanowarp

import "gonum.org/v1/gonum/dsp/fourier"

// Masuyama, Yoshiki, Kohei Yatabe, and Yasuhiro Oikawa.
// "Griffin–Lim like phase recovery via alternating direction method of multipliers."
// IEEE Signal Processing Letters 26.1 (2018): 184-188.
// https://ieeexplore.ieee.org/iel7/97/4358004/08552369.pdf

type admm struct {
	fft    *fourier.FFT
	s      []float64
	norm   float64
	hop    int
	Wf, Wr []float64
}

type stftTensor struct {
	data    []complex128
	shape   [2]int
	padding int
	view    [][]complex128
}

func (n *admm) stft(dst [][]complex128, sa []float64) [][]complex128 {
	nfft := n.fft.Len()
	for i := 0; i < len(sa); i += n.hop {
		copy(n.s[max(0, nfft/2-i):], sa[max(0, i-nfft/2):min(len(sa), i+nfft/2)])
		mul(n.s, n.Wf)
		n.fft.Coefficients(dst[i/n.hop], n.s)
	}
	return dst
}

func (n *admm) istft(dst []float64, fr [][]complex128) []float64 {
	nfft := n.fft.Len()
	for j := range fr {
		i := j * n.hop
		n.fft.Sequence(n.s, fr[j])
		mul(n.s, n.Wr)
		scale(n.s, 1/n.norm)
		add(dst[max(0, i-nfft/2):min(len(dst), i+nfft/2)], n.s)
	}
	return dst
}

func (n *admm) project(dst [][]complex128) [][]complex128 {
	return n.stft(dst, n.istft(nil, dst))
}
