package nanowarp

import (
	"fmt"
	"os"

	"gonum.org/v1/gonum/cmplxs"
	"gonum.org/v1/gonum/floats"
)

func (n *warper) processGla(lin, rin, lout, rout []float64, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).processGla`)
	// println := func(a ...any) {}

	input := make2[float64](2, len(lin))
	grainbuf := make2[float64](2, n.nfft)
	ingrains := make3[float64](n.lah, 2, n.nfft)
	outgrains := make3[float64](n.lah, 2, n.nfft)
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

	a := &n.a
	for j := -n.nbuf; ; j += n.hop * (n.lah - 1) {
		if j > len(lout)-n.nbuf {
			break
		}
		l := min((len(lout)-j)/n.hop, n.lah)
		for i := range l {
			j := j + i*n.hop
			getgrain(ingrains[i], j)
		}

		for t := range ingrains[:l] {
			n.analyze2(ingrains[:l][t], a.Cs[t+n.lah], a.Xs[t+n.lah], a.Ms[t+n.lah])
		}

		for i := range outgrains {
			n.synthesize2(outgrains[i], a.Cs[i], a.Xs[i])
		}

		g := outgrains[1:l]
		for i := range g {
			addgrain(j+(i+1)*n.hop, g[i])
		}
	}
}

func (n *warper) analyze2(present [][]float64, C [][]complex128, X []complex128, M []float64) {
	a := &n.a

	clear(a.Mid)
	for ch := range present {
		floats.Add(a.Mid, present[ch])
		n.enfft(C[ch], a.W, present[ch])
	}

	n.enfft(X, a.W, a.Mid)

	for w := range a.X {
		M[w] = mag(X[w])
	}
	for ch := range present {
		for w := range a.X {
			C[ch][w] /= X[w]
		}
	}
}

func (n *warper) synthesize2(output [][]float64, C [][]complex128, P []complex128) {
	a := &n.a
	for ch := range output {
		cmplxs.MulTo(a.X, C[ch], P)
		n.defft(output[ch], a.X, true)
	}
}
