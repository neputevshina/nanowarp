package nanowarp

import (
	"math"
	"slices"

	"gonum.org/v1/gonum/dsp/fourier"
)

// hann is a Hann window function.
//
// Dual to itself (achieves COLA/perfect reconstruction) at any integer overlap greater than 2.
func hann(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out)-1)
		out[i] = 0.5 * (1 - math.Cos(2*math.Pi*x))
	}
}

// blackmanHarrisClassic is a, counterintuitively, three-term Nutall minimum sidelobe window function
// also known as “three-term Blackman-Harris window”[1].
//
// Dual to itself (achieves COLA/perfect reconstruction) at overlap 6.
//
// Fun fact: frequency response in the article “Window functon” in English Wikipedia is for
// four-term -92 dB Blackman-Harris. And the picture of it was created and uploaded by Olli Niemitalo,
// the creator of asymmetric window function defined below. Correct parameters for the pictured window
// can be found in German version of this article. Still sounds better than four-term version.
//
// [1]: Doerry, Armin W. Catalog of window taper functions for sidelobe control. No. SAND2017-4042.
// Sandia National Laboratories (SNL-NM), Albuquerque, NM (United States), 2017.
// https://www.osti.gov/servlets/purl/1365510
func blackmanHarrisClassic(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out)-1)
		out[i] = .4243801 - .4973406*math.Cos(2*math.Pi*x) + .0782793*math.Cos(4*math.Pi*x)
	}
}

// blackmanHarris92dB is a four-term -92 dB Blackman-Harris window function.
// See blackmanHarrisClassic.
func blackmanHarris92dB(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out)-1)
		out[i] = .35875 - .48829*math.Cos(2*math.Pi*x) + .14128*math.Cos(4*math.Pi*x) - .01168*math.Cos(6*math.Pi*x)
	}
}

// avciNacaroglu is an Avci-Nacaroglu window function [1].
// It is related to DPSS window like Kaiser, but doesn't require
// modified Bessel function to calculate.
// Duality to itself is achieved by tuning the parameter a.
// For example, COLA at overlap 4 can be approximately achieved with a = 1.78.
//
// [1]: Doerry, Armin W. Catalog of window taper functions for sidelobe control. No. SAND2017-4042.
// Sandia National Laboratories (SNL-NM), Albuquerque, NM (United States), 2017.
// https://www.osti.gov/servlets/purl/1365510
func avciNacaroglu(out []float64, a float64) {
	for i := range out {
		x := float64(i) / float64(len(out)-1)
		out[i] = math.Exp(math.Pi*a*math.Sqrt(1-(4*(x-0.5)*(x-0.5))) - 1)
	}
	m := slices.Max(out)
	for i := range out {
		out[i] /= m
	}
}

// niemitalo is an asymmetric windowing function, designed by Olli Niemitalo.
//
// See https://dsp.stackexchange.com/questions/2337/fft-with-asymmetric-windowing
func niemitalo(out []float64) {
	nfft := float64(len(out))
	clear(out)
	sin, cos := math.Sin, math.Cos
	for i := nfft / 4; i < nfft*7/8; i++ {
		x := 2 * math.Pi * ((i+0.5)/nfft - 1.75)
		out[int(i)] = 2.57392230162633461887 - 1.58661480271141974718*cos(x) + 3.80257516644523141380*sin(x) -
			1.93437090055110760822*cos(2*x) - 3.27163999159752183488*sin(2*x) + 3.26617449847621266201*cos(3*x) -
			0.30335261753524439543*sin(3*x) - 0.92126091064427817479*cos(4*x) + 2.33100177294084742741*sin(4*x) -
			1.19953922321306438725*cos(5*x) - 1.25098147932225423062*sin(5*x) + 0.99132076607048635886*cos(6*x) -
			0.34506787787355830410*sin(6*x) - 0.04028033685700077582*cos(7*x) + 0.55461815542612269425*sin(7*x) -
			0.21882110175036428856*cos(8*x) - 0.10756484378756643594*sin(8*x) + 0.06025986430527170007*cos(9*x) -
			0.05777077835678736534*sin(9*x) + 0.00920984524892982936*cos(10*x) + 0.01501989089735343216*sin(10*x)
	}
	for i := 0; i < int(nfft)/8; i++ {
		nfft := int(nfft)
		out[nfft-1-i] = (1 - out[nfft*3/4-1-i]*out[nfft*3/4+i]) / out[nfft/2+i]
	}
	copy(out, out[int(nfft)*2/8:])
	clear(out[int(nfft)*6/8:])
}

// windowGain returns the squared window gain for correcting the output grain
// gain after double (for fft and ifft) application of the window.
func windowGain(w []float64) (a float64) {
	for _, e := range w {
		a += e * e
	}
	a /= float64(len(w))
	return
}

// windowT calculates the time ramp multiplied window
// function for time-frequency reassignment.
func windowT(out, w []float64) {
	n := float64(len(w))
	for i := range w {
		out[i] = w[i] * mix(-n/2, n/2+1, float64(i)/n)
	}
}

// windowDx calculates the derivative of w using DFT.
func windowDx(out, w []float64) {
	f := fourier.NewFFT(len(w))
	s := f.Coefficients(nil, w)
	for i := range s {
		s[i] *= complex(0, float64(i))
		s[i] /= complex(float64(len(w)), 0)
	}
	s[len(w)/2] = 0
	f.Sequence(out, s)
	for i := range out {
		out[i] = -out[i] * math.Pi * 2
	}
}

// windowDualUniform calculates dual window for uniform synthesis hop size
// from a given symmetrical window to satisfy constant overlap-add (COLA)
// condition for perfect reconstruction.
//
// It returns the total gain after application of both windows.
//
// It can't necessarily be used to generate windows for NSDGT.
func windowDualUniform(out, w []float64, hop int) (gain float64) {
	buf := make([]float64, 4*len(w))
	w2 := slices.Clone(w)
	mul(w2, w)
	half := len(w) / 2
	for i := half; i < len(buf)-half; i += hop {
		add(buf[i-half:i+half], w2)
	}
	copy(out, w)
	mul(out, buf[len(w):2*len(w)])
	return sum(mul(slices.Clone(w), out)) / float64(len(w))
}
