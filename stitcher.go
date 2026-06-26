package nanowarp

import (
	"fmt"
	"math/cmplx"
	"os"

	"gonum.org/v1/gonum/cmplxs"
	"gonum.org/v1/gonum/floats"
)

func (n *warper) processGla(lin, rin, lout, rout []float64, phasor []float64, stretch float64) {
	fmt.Fprintln(os.Stderr, `(*warper).processGla`)
	// println := func(a ...any) {}

	lah := 2 * n.olap

	input := make2[float64](2, len(lin))
	grainbuf := make2[float64](2, n.nfft)
	ingrains := make3[float64](lah, 2, n.nfft)
	outgrain := make2[float64](2, n.nfft)
	copy(input[0], lin)
	copy(input[1], rin)

	getgrain := func(ingrain [][]float64, j int) {
		i := int(phasor[max(0, j)])
		if i > len(lin)-n.nbuf {
			return
		}
		for ch := range grainbuf {
			// fill(ingrain[ch], 1)
			// continue
			clear(ingrain[ch])
			copy(ingrain[ch][max(0, -i+n.nbuf/2):], input[ch][clamp(0, len(lout)-n.nfft, i-n.nbuf/2):])
		}
	}
	addgrain := func(j int, grainbuf [][]float64) {
		loutgrain := lout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		add(loutgrain, grainbuf[0][clamp(0, n.nbuf, -j):])

		routgrain := rout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		add(routgrain, grainbuf[1][clamp(0, n.nbuf, -j):])
	}

	known := make([]bool, lah*2)
	fill(known[:lah], true)
	known[lah*2-1] = true
	// .  .  .  .  ?  ?  ?  .
	// [________]  [_____]  _
	// lah         lah-1    1

	a := &n.a
	for j := -n.nbuf; ; j += n.hop * lah {
		// println(float64(j) / float64(len(lout)))
		if j > len(lout)-n.nbuf {
			break
		}
		l := min((len(lout)-j)/n.hop, lah)
		for t := range l {
			j := j + t*n.hop
			getgrain(ingrains[t], j)
		}

		for t := range ingrains[:l] {
			n.analyze2(ingrains[:l][t], a.Cs[t], a.Xs[lah+t], a.Phs[t], a.Ms[lah+t], a.Fadvs[t+1], a.Tadvs[t+1], stretch)
		}

		// n.integrate(a.Fadvs, a.Tadvs, a.Ms[lah-1:], a.Phs[:lah+1], n.arm, known[lah-1:])
		n.admm.admm(a.Xs[:2*lah], a.Ms[:2*lah], known, 15, 0.1)

		for t := range ingrains[:l] {
			n.synthesize2(outgrain, a.Cs[t], a.Xs[lah+t])
			apply2(outgrain, bitsafe)

			addgrain(j+t*n.hop, outgrain)
			copy(a.Xs[t], a.Xs[lah+t])
			copy(a.Ms[t], a.Ms[lah+t])
		}
		// waveform.Dump(nil, a.sbig[:a.hop*8+2*a.nbuf])

		// println()
		// println()
		// println()
		// println()
	}
}

func (n *warper) analyze2(present [][]float64, C [][]complex128, X []complex128, Ph, M, Fadv, Tadv []float64, speedup float64) {
	a := &n.a

	clear(a.Mid)
	for ch := range present {
		floats.Add(a.Mid, present[ch])
		clear(a.S)
		copy(a.S, present[ch])
		mul(a.S, a.W)
		n.fft.Coefficients(C[ch], a.S)
	}

	copy(a.S, a.Mid)
	mul(a.S, a.W)
	n.fft.Coefficients(X, a.S)
	n.enfft(a.Xd, a.Wd, a.Mid)
	n.enfft(a.Xt, a.Wt, a.Mid)

	for w := range a.X {
		Fadv[w] = princarg(getfadv(X, a.Xt, 2./n.osamp/speedup)(w))
		Tadv[w] = gettadv(X, a.Xd, float64(n.olap)*n.osamp)(w)

	}

	for w := range a.X {
		M[w], Ph[w] = cmplx.Polar(X[w])
	}
	for ch := range present {
		for w := range a.X {
			C[ch][w] /= X[w]
		}
	}
}

func (n *warper) synthesize2(output [][]float64, C [][]complex128, X []complex128) {
	a := &n.a
	for ch := range output {
		cmplxs.MulTo(a.X, C[ch], X)
		n.fft.Sequence(a.S, X)
		mul(a.S, a.W)
		floats.Scale(1./n.norm, a.S)
		copy(output[ch], a.S)
	}
}
