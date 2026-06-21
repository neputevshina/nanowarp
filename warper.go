package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"slices"

	"github.com/neputevshina/nanowarp/oscope"
	"gonum.org/v1/gonum/cmplxs"
	"gonum.org/v1/gonum/dsp/fourier"
	"gonum.org/v1/gonum/floats"
)

type warper struct {
	nfft  int     // DFT size, a power of 2
	nbuf  int     // Effective window size, nbuf<nfft
	hop   int     // Window output hop size
	lah   int     // Non-causal PGHI lookahead in frames (hop sizes)
	nbins int     // nfft/2+1, Number of DFT bins
	olap  int     // nbuf/hop, Window ovelap
	osamp float64 // nfft/nbuf, Zero-padding ratio

	root *Nanowarp

	fft         *fourier.FFT
	arm         [][]bool // PGHI done mask
	norm, wgain float64  // Global normalization factor and grain-only normalization factor
	heap        hp       // PGHI heap

	a wbufs
}

type wbufs struct {
	S, Mid, M, P, F            []float64 // Scratch buffers
	Ph                         []float64 `size:"nbins"` // Current phase
	Past, Future               []float64 // Phase accumulators
	Fadv, Tadv                 []float64
	W, Wr, Wd, Wt, Wdt         []float64      // Window functions
	X, Y, Xd, Xt, L, R, Lo, Ro []complex128   // Complex spectra
	C, Co                      [][]complex128 // Channels

	Cs                        [][][]complex128
	Phs, Fadvs, Tadvs, Ms, As [][]float64 `size:"lah"`
}

func warperNew(nbuf, osamp, olap, nch int, nanowarp *Nanowarp) (n *warper) {
	// FIXME Only 2x oversampling works, no more, no less.
	nfft := nextpow2(nbuf * osamp)
	n = &warper{
		nfft:  nfft,
		nbins: nfft/2 + 1,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		olap:  olap,
		osamp: float64(osamp),
		root:  nanowarp,
		lah:   olap*300 + 1,
	}
	a := &n.a

	makeslices(a, n.nbins, nfft, nch, n.lah)

	n.arm = make([][]bool, n.lah)
	for i := range n.arm {
		n.arm[i] = make([]bool, n.nbins)
	}

	s := func(w []float64) []float64 {
		// return w[nfft/2-nbuf/2 : nfft/2+nbuf/2]
		return w[:nbuf]
	}
	blackmanHarris(s(a.W))

	// FIXME Destination is the first argument by convention.
	windowDx(s(a.W), s(a.Wd))
	windowT(s(a.W), s(a.Wt))
	windowT(s(a.Wd), s(a.Wdt))

	copy(s(a.Wr), s(a.W))
	slices.Reverse(s(a.Wr))
	n.wgain = windowGain(n.a.W)
	n.norm = float64(nfft) * float64(n.olap) * n.osamp * n.wgain

	// waveform.Dump(nil, a.W)
	// waveform.Dump(nil, a.Wt)
	// waveform.Dump(nil, a.Wd)

	n.fft = fourier.NewFFT(nfft)
	n.heap = make(hp, n.lah*n.nbins) // 2 for future and past.

	return
}

func (n *warper) process3old(lin, rin, lout, rout []float64, coeffs, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	get := func() []float64 { return make([]float64, n.nfft) }
	lgrainbuf, rgrainbuf := get(), get()
	lingrain, ringrain := get(), get()
	lfuture, rfuture := get(), get()

	lastone := 0
	for j := -n.nbuf; ; j += n.hop {
		if j > len(lout)-n.nbuf {
			break
		}
		fi := func(j int) int { return int(phasor[max(0, j)]) }
		i := fi(j)
		c := 1 / coeffs[max(0, j)]

		cut := func(lingrain, lin []float64, i int) {
			copy(lingrain[max(0, -(i-n.nbuf/2)):], lin[max(0, (i-n.nbuf/2)):clamp(0, len(lin), i+n.nbuf/2)])
		}
		cut(lingrain, lin, i)
		cut(ringrain, rin, i)
		if c == 1 {
			cut(lfuture, lin, fi(j+n.hop))
			cut(rfuture, rin, fi(j+n.hop))
			add(lfuture, rfuture)
			scale(lfuture, 0.5)
		}

		n.advanceOld(lingrain, ringrain, lgrainbuf, rgrainbuf, lfuture, abs(c), c == 1)

		// Cut pre-echo in transient regions.
		d := j - lastone
		if c == 1 {
			lastone = j
		} else if d < n.nbuf/2 {
			z := func(grain []float64) {
				rr := grain[max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grain[:n.nbuf/2-d-n.hop], 0)
			}
			z(lgrainbuf)
			z(rgrainbuf)
		}

		if n.root.opts.Onsets && c != 1 {
			clear(lgrainbuf)
			clear(rgrainbuf)
		}

		loutgrain := lout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		routgrain := rout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]

		add(loutgrain, lgrainbuf[clamp(0, n.nbuf, -j):])
		add(routgrain, rgrainbuf[clamp(0, n.nbuf, -j):])
	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advanceOld(lingrain, ringrain, loutgrain, routgrain, future []float64, stretch float64, reset bool) {
	a := &n.a
	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	for i := range lingrain {
		a.Mid[i] = lingrain[i] + ringrain[i]
	}

	enfft(a.L, a.W, lingrain)
	enfft(a.R, a.W, ringrain)

	enfft(a.X, a.W, a.Mid)
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	// TODO Works ONLY with 2x FFT oversampling.
	olap := float64(n.nbuf / n.hop)
	osampc := 2.
	tadv := gettadv(a.X[:n.nbins], a.Xd, osampc*olap)
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch)

	for w := range a.X {
		// TODO Probably it will be more numerically stable to limit the phase accuum to
		// 0..1 and scale back to -π..π range at the poltocar conversion.
		// fadv must return 0..1 accordingly, simply defer π multiplication till the end.
		a.Fadv[w] = princarg(fadv(w))
	}

	for w := range a.X {
		// Encode stereo phase differences and stretch mid only, keep original magnitudes.
		// NB: Phase difference in polar coordinates is complex division in cartesian.
		//     Phase sum is conversely a multiply.
		//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
		//
		// See “Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
		// real-time realization for a commercial music product.”
		m := mag(a.X[w])
		p := a.X[w] / complex(m, 0)
		if m < 1e-6 {
			p = complex(1, 0)
		}
		a.M[w] = m

		a.L[w] /= p
		a.R[w] /= p
	}

	n.heap = make(hp, n.nbins)

	clear(a.F)
	if reset {
		enfft(a.Y, a.W, future)
	}
	clear(n.arm[0])
	for w := range a.X {
		n.arm[0][w] = true
		n.heap[w] = heaptriple{a.P[w], w, -1}
	}

	heapInit(&n.heap)

	// Do PGHI.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[0][w] {
				a.F[w]++
				adv := tadv(w)
				a.Ph[w] = a.Past[w] + adv
				n.arm[0][w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, 0})
			}
		case 0:
			if w > 1 && n.arm[0][w-1] {
				a.F[w]++
				adv := -a.Fadv[w-1]
				a.Ph[w-1] = a.Ph[w] + adv
				n.arm[0][w-1] = false
				heapPush(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[0][w+1] {
				a.F[w]++
				adv := a.Fadv[w+1]
				a.Ph[w+1] = a.Ph[w] + adv
				n.arm[0][w+1] = false
				heapPush(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
			}
		}
	}

	copy(a.P, a.M)
	for w := range a.Ph {
		// Add stereo phase differences back through multiplication.
		ph := a.Ph[w]
		if reset && a.F[w] < 3 {
			ph = cmplx.Phase(a.X[w])
		}
		phasor := cmplx.Rect(1, ph)
		a.Lo[w] = a.L[w] * phasor
		a.Ro[w] = a.R[w] * phasor
		a.Past[w] = princarg(ph)
	}
	goto skip
skip:
	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.Wr)
		copy(out, a.S)
	}
	defft(loutgrain, a.Lo)
	defft(routgrain, a.Ro)
}

func (n *warper) process3(lin, rin, lout, rout []float64, coeffs, phasor []float64, causal bool) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	println := func(a ...any) {}

	input := make2(2, len(lin))
	grainbuf := make2(2, n.nfft)
	ingrain := make2(2, n.nfft)
	copy(input[0], lin)
	copy(input[1], rin)
	lastone, firstone := 0, 0
	latch := true
	c := func(j int) float64 {
		c := 1.
		if j > 0 {
			c = coeffs[j]
			if c != c || math.IsInf(c, 0) || c == 0 {
				c = 1
			}
		}
		return c
	}

	for j := -n.nbuf; ; j += n.hop {
		if j > len(lout)-n.nbuf {
			break
		}

		if !causal && j > 0 {
			// TODO Compensate this.
			future := j + n.hop*n.lah
			filt := future-firstone > n.root.fs*n.root.opts.TransientMs/1000
			if future < len(lout) && coeffs[future] == 1 && filt {
			}
		}

		d := j - lastone
		if c(j) == 1 {
			lastone = j
			if latch {
				println(`firstone:`, j)
				firstone = j
			}
			latch = false
		} else {
			latch = true
		}

		i := int(phasor[max(0, j)])

		if i > len(lin)+n.nbuf/2 {
			break
		}
		for ch := range grainbuf {
			clear(ingrain[ch])
			copy(ingrain[ch][max(0, -i+n.nbuf/2):], input[ch][max(0, i-n.nbuf/2):])
		}

		if n.root.opts.Quality >= 0 && c(j) == 1 {
			println(`reset past:`, j)
			n.resetPast(ingrain)
			n.bypassGrain(ingrain, grainbuf)
		} else {
			n.advance([][][]float64{ingrain}, grainbuf, c(j))
		}

		// Cut pre-echo in transient regions.
		// TODO Works strange after implementing non-causality.
		if c(j) != 1 && d < n.nbuf/2 {
			for ch := range grainbuf {
				rr := grainbuf[ch][max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grainbuf[ch][:n.nbuf/2-d-n.hop], 0)
			}
		}

		if n.root.opts.Onsets && c(j) != 1 {
			for ch := range grainbuf {
				clear(grainbuf[ch])
			}
		}

		loutgrain := lout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		add(loutgrain, grainbuf[0][clamp(0, n.nbuf, -j):])

		routgrain := rout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		add(routgrain, grainbuf[1][clamp(0, n.nbuf, -j):])
	}
}

func (n *warper) process5(lin, rin, lout, rout []float64, coeffs, phasor []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process4`)
	// println := func(a ...any) {}

	input := make2(2, len(lin))
	grainbuf := make2(2, n.nfft)
	ingrains := make3(n.lah, 2, n.nfft)
	outgrains := make3(n.lah, 2, n.nfft)
	copy(input[0], lin)
	copy(input[1], rin)
	speedups := make([]float64, n.lah)

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
	_ = addgrain

	for j := -n.nbuf; ; j += n.hop * (n.lah - 1) {
		if j > len(lout)-n.nbuf {
			break
		}
		l := min((len(lout)-j)/n.hop, n.lah)
		println(l)
		for i := range l {
			j := j + i*n.hop
			getgrain(ingrains[i], j)
			speedups[i] = coeffs[max(0, j)]
		}
		n.stitch(true, true, ingrains[:l], outgrains, speedups)
		g := outgrains[1:l]
		for i := range g {
			addgrain(j+(i+1)*n.hop, g[i])
		}
	}
}

func (n *warper) enfft(x []complex128, w, grain []float64) {
	a := &n.a
	clear(a.S)
	copy(a.S, grain)
	if w != nil {
		mul(a.S, w)
	}
	n.fft.Coefficients(x, a.S)
}

func (n *warper) defft(out []float64, x []complex128, w bool) {
	a := &n.a
	n.fft.Sequence(a.S, x)
	floats.Scale(1./n.norm, a.S)
	if w {
		mul(a.S, a.Wr)
	}
	copy(out, a.S)
}

func (n *warper) resetPast(present [][]float64) {
	a := &n.a
	clear(a.Mid)
	for ch := range present {
		floats.Add(a.Mid, present[ch])
		n.enfft(a.C[ch], a.W, present[ch])
	}
	n.enfft(a.X, a.W, a.Mid)
	for w := range a.Ph {
		a.P[w], a.Past[w] = cmplx.Polar(a.X[w])
	}
}

func (n *warper) bypassGrain(present, output [][]float64) {
	a := &n.a
	for ch := range present {
		// This is the only known way to correctly scale gains.
		n.enfft(a.Co[ch], a.W, present[ch])
		n.defft(output[ch], a.Co[ch], true)
		// And this works with -40 dB difference.
		// copy(output[ch], present[ch])
		// mul(output[ch], a.W)
		// mul(output[ch], a.Wr)
		// scale(output[ch], n.wgain*float64(n.nbuf)/float64(n.hop))
		// Probably there is something wrong in gonum FFT implementation.
	}
}

func (n *warper) bash(present [][]float64, C [][]complex128, M, Ph []float64, Mid []float64) {
	a := &n.a
	clear(Mid)
	for ch := range present {
		floats.Add(Mid, present[ch])
		n.enfft(C[ch], a.W, present[ch])
	}
	n.enfft(a.X, a.W, Mid)
	for w := range n.nbins {
		M[w], Ph[w] = cmplx.Polar(a.X[w])
	}
}

func (n *warper) analyze(present [][]float64, C [][]complex128, Fadv, Tadv, M, Mid []float64, speedup float64) {
	a := &n.a

	clear(Mid)
	for ch := range present {
		floats.Add(Mid, present[ch])
		n.enfft(C[ch], a.W, present[ch])
	}

	n.enfft(a.X, a.W, Mid)
	n.enfft(a.Xd, a.Wd, Mid)
	n.enfft(a.Xt, a.Wt, Mid)
	n.enfft(a.Y, a.Wdt, Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	for w := range a.X {
		// TODO Probably it will be more numerically stable to limit the phase accuum to
		// 0..1 and scale back to -π..π range at the poltocar conversion.
		// fadv must return 0..1 accordingly, simply defer π multiplication till the end.
		Fadv[w] = princarg(getfadv(a.X, a.Xt, 2./n.osamp/speedup)(w))
		Tadv[w] = gettadv(a.X, a.Xd, float64(n.olap)*n.osamp)(w)
	}

	for w := range a.X {
		// Encode stereo phase differences and stretch mid only, keep original magnitudes.
		// NB: Phase difference in polar coordinates is complex division in cartesian.
		//     Phase sum is conversely a multiply.
		//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
		//
		// See “Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
		// real-time realization for a commercial music product.”
		m := mag(a.X[w])
		p := a.X[w] / complex(m, 0)
		if m < 1e-6 {
			p = complex(1, 0)
		}
		// 4.5 dB/octave tilt
		M[w] = m //* (1 - math.Pow(float64(w/n.nbins), 1./2.82))
		a.Y[w] = p
	}
	for ch := range present {
		for w := range a.X {
			C[ch][w] /= a.Y[w]
		}
	}
}

func (n *warper) integrate(Fadv, Tadv, M [][]float64, Ph [][]float64, arm [][]bool) {
	n.heap = n.heap[:0]
	clear(n.heap)

	// Frames on sides of the framebuffer are the past frame after the
	// previous transient, and the future frame, which is the first
	// frame of a transient.
	//
	// They are ground truths for the current step of integration, so they are
	// added to the heap and not armed for phase reconstruction.
	for t := range Ph {
		if t > 0 && t < n.lah {
			fill(arm[t], true)
		} else {
			fill(arm[t], false)
			for w := range n.nbins {
				n.heap = append(n.heap, heaptriple{M[t][w], w, t})
			}
		}
	}
	heapInit(&n.heap)

	oscope.Enable = true

	var rights, ups, downs, lefts [][]float64
	if oscope.Enable {
		rights = make2(len(Fadv), n.nbins)
		ups = make2(len(Fadv), n.nbins)
		downs = make2(len(Fadv), n.nbins)
		lefts = make2(len(Fadv), n.nbins)
	}

	// Do PGHI with causal von Neumann neighborhood.
	hood := [][]float64{
		{0, 1, 0},
		{0, 0, 1},
		{0, 1, 0}}
	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		for w := range hood {
			for t, v := range hood[w] {
				type f = float64
				tn, wn := t-(len(hood)-1)/2, w-(len(hood[t])-1)/2
				t, w := h.t+tn, h.w+wn
				// Check bounds.
				if clamp(0, len(arm)-1, t) != t || clamp(0, len(arm[t])-1, w) != w {
					continue
				}
				// Skip complete.
				if !arm[t][w] || v == 0 {
					continue
				}
				Ph[t][w] = Ph[h.t][h.w] + f(tn)*Tadv[t][w] + f(wn)*Fadv[t][w]
				arm[t][w] = false
				heapPush(&n.heap, heaptriple{M[t][w], w, t})
			}
		}
	}

	// for len(n.heap) > 0 {
	// 	h := heapPop(&n.heap)
	// 	w, t := h.w, h.t
	// 	// if t >= 1 && arm[t-1][w] {
	// 	// 	if oscope.Enable {
	// 	// 		lefts[t][w] = 1
	// 	// 	}
	// 	// 	Ph[t-1][w] = Ph[t][w] - Tadv[t-1][w]
	// 	// 	arm[t-1][w] = false
	// 	// 	heapPush(&n.heap, heaptriple{M[t-1][w], w, t - 1})
	// 	// }
	// 	if t < len(Ph)-1 && arm[t+1][w] {
	// 		if oscope.Enable {
	// 			// rights[t+1][w] = 1
	// 			rights[t][w] = 1
	// 		}
	// 		Ph[t+1][w] = Ph[t][w] + Tadv[t+1][w]
	// 		arm[t+1][w] = false
	// 		heapPush(&n.heap, heaptriple{M[t+1][w], w, t + 1})
	// 	}
	// 	if w >= 1 && arm[t][w-1] {
	// 		if oscope.Enable {
	// 			downs[t][w] = 1
	// 			// downs[t][w-1] = 1
	// 		}
	// 		Ph[t][w-1] = Ph[t][w] - Fadv[t][w-1]
	// 		arm[t][w-1] = false
	// 		heapPush(&n.heap, heaptriple{M[t][w-1], w - 1, t})
	// 	}
	// 	if w < n.nbins-1 && arm[t][w+1] {
	// 		if oscope.Enable {
	// 			ups[t][w] = 1
	// 			// ups[t][w+1] = 1
	// 		}
	// 		Ph[t][w+1] = Ph[t][w] + Fadv[t][w+1]
	// 		arm[t][w+1] = false
	// 		heapPush(&n.heap, heaptriple{M[t][w+1], w + 1, t})
	// 	}
	// }
	if oscope.Enable {
		for _, e := range lefts {
			oscope.Oscope(slices.Clone(e), oscope.Name(`lefts`))
		}
		for _, e := range rights {
			oscope.Oscope(slices.Clone(e), oscope.Name(`rights`))
		}
		for _, e := range ups {
			oscope.Oscope(slices.Clone(e), oscope.Name(`ups`))
		}
		for _, e := range downs {
			oscope.Oscope(slices.Clone(e), oscope.Name(`downs`))
		}
		// for t := range downs {
		// 	for w := range downs[t] {
		// 		n.a.As[t][w] = boolfloat(rights[t][w]+ups[t][w]+downs[t][w] <= 1)
		// 	}
		// }

		// apply2(M, bitsafe)
		// apply2(M, atodb)
		// apply2(M, bitsafe)
		for _, e := range M {
			oscope.Oscope(slices.Clone(e), oscope.Name(`mag`))
		}
		for t, e := range M {
			clear(n.a.S)
			for w := 1; w < len(e)-1; w++ {
				// ups[t][w] = math.Log10(bitsafe(e[w] / e[w+1]))
				// downs[t][w] = -ups[t][w]
				// rights[t][w] = math.Log10(bitsafe(e[w] / M[min(len(M)-1, t+1)][w]))

				if e[w-1] < e[w] {
					downs[t][w] = 1
				}
				if e[w+1] < e[w] {
					ups[t][w] = 1
				}
				if M[min(len(M)-1, t+1)][w] < e[w] {
					rights[t][w] = 1
				}
			}
			oscope.Oscope(slices.Clone(n.a.S), oscope.Name(`mags`))
			clear(n.a.S)
			for w := 1; w < len(e)-1; w++ {
				if downs[t][w] == 1 && ups[t][w] == 1 {
					n.a.S[w] = 1
				}
			}
			oscope.Oscope(slices.Clone(n.a.S), oscope.Name(`ridge2`))
		}
		histonorm2(rights)
		histonorm2(ups)
		histonorm2(downs)
		for _, e := range rights {
			clear(e)
			oscope.Oscope(slices.Clone(e), oscope.Name(`rights-mag`))
		}
		for _, e := range ups {
			clear(e)
			oscope.Oscope(slices.Clone(e), oscope.Name(`ups-mag`))
		}
		for _, e := range downs {
			clear(e)
			oscope.Oscope(slices.Clone(e), oscope.Name(`downs-mag`))
		}
		for t := range downs {
			n.partition(n.a.As[t], M[t], 5)
			oscope.Oscope(slices.Clone(n.a.As[t][:n.nbins]), oscope.Name(`ridges-mag`))
		}
	}
}

func (n *warper) partition(output []float64, input []float64, c int) {
	if len(input) < c {
		return
	}
	i := argmax(input)
	output[i] = 1
	if i > c+1 && i < len(input)-c-3 {
		n.partition(output[:i-c], input[:i-c], c-1)
		n.partition(output[i+c:], input[i+c:], c+1)
	}
}

func (n *warper) synthesize(output [][]float64, C [][]complex128, Ph, Mask []float64) {
	a := &n.a
	for w := range Ph {
		// Add stereo phase differences back through complex multiplication.
		m := 1.
		if Mask != nil {
			m = Mask[w]
		}
		a.Y[w] = cmplx.Rect(m, Ph[w])
	}
	for ch := range output {
		cmplxs.MulTo(a.X, C[ch], a.Y)
		n.defft(output[ch], a.X, true)
	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(fs [][][]float64, output [][]float64, speedup float64) {
	a := &n.a

	present := fs[0]

	n.analyze(present, a.C, a.Fadv, a.Tadv, a.M, a.Mid, speedup)

	n.integrate([][]float64{nil, a.Fadv}, [][]float64{nil, a.Tadv}, [][]float64{a.P, a.M}, [][]float64{a.Past, a.Ph}, n.arm)

	for w := range a.Ph {
		a.Past[w] = princarg(a.Ph[w])
	}
	copy(a.P, a.M)

	n.synthesize(output, a.C, a.Ph, nil)
}

// stitch reconstructs the phase between fs[0] and fs[n.lah-1] frames.
func (n *warper) stitch(g0, glast bool, fs [][][]float64, output [][][]float64, speedup []float64) {
	a := &n.a

	last := len(fs) - 1
	if g0 {
		n.bash(fs[0], a.Cs[0], a.Ms[0], a.Phs[0], a.Mid)
	}
	if glast {
		n.bash(fs[last], a.Cs[last], a.Ms[last], a.Phs[last], a.Mid)
	}
	for t := range last - 1 {
		n.analyze(fs[t+1], a.Cs[t+1], a.Fadvs[t+1], a.Tadvs[t+1], a.Ms[t+1], a.Mid, speedup[t+1])
		// println(t + 1)
		// waveform.Dump(nil, fs[t+1][0])
	}

	n.integrate(a.Fadvs, a.Tadvs, a.Ms, a.Phs, n.arm)

	for i := range last + 1 {
		n.synthesize(output[i], a.Cs[i], a.Phs[i], nil)
	}
}

func getfadv(x, xt []complex128, stretch float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) == 0 {
			return 0
		}
		// NOTE Try len(x)-1 instead. Sounds worse on my $4 speakers.
		// FIXME Works ONLY with nbuf=4096, nfft=8192 (oversampling 2).
		return -real(xt[j]/x[j])/float64(len(x))*math.Pi*stretch - math.Pi/2
	}
}

func gettadv(x, xd []complex128, olap float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) < 1e-6 {
			return 0
		}
		return (math.Pi*float64(j) + imag(xd[j]/x[j])) / (olap / 2)
	}
}
