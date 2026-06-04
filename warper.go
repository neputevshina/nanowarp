package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"slices"

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
	arm         []bool  // PGHI done mask, 2D, shift addressed
	norm, wgain float64 // Global normalization factor and grain-only normalization factor
	heap        hp      // PGHI heap

	a wbufs
}

type wbufs struct {
	S, Mid, M, P, F            []float64 // Scratch buffers
	Ph                         []float64 `size:"nbins"` // Current phase
	Past, Future               []float64 // Phase accumulators
	Fadv                       []float64
	W, Wr, Wd, Wt              []float64      // Window functions
	X, E, Xd, Xt, L, R, Lo, Ro []complex128   // Complex spectra
	C, Co                      [][]complex128 // Channels
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
		lah:   olap * osamp * 2,
	}
	a := &n.a

	makeslices(a, n.nbins, nfft, nch)
	n.arm = make([]bool, n.nbins*n.lah)

	s := func(w []float64) []float64 {
		// return w[nfft/2-nbuf/2 : nfft/2+nbuf/2]
		return w[:nbuf]
	}
	blackmanHarris(s(a.W))

	windowDx(s(a.W), s(a.Wd))
	windowT(s(a.W), s(a.Wt))
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

func (n *warper) process3(lin, rin, lout, rout []float64, coeffs, phasor []float64, causal bool) {
	fmt.Fprintln(os.Stderr, `(*warper).process3`)
	println := func(a ...any) {}

	input := make2(2, len(lin))
	grainbuf := make2(2, n.nfft)
	ingrain := make2(2, n.nfft)
	copy(input[0], lin)
	copy(input[1], rin)
	h := n.hop
	pin, lastone, firstone := 0, 0, 0
	latch := true
	for j := -n.nbuf; ; j += h {
		if j > len(lout)-n.nbuf {
			break
		}

		c := func() float64 {
			c := 1.
			if j > 0 {
				c = 1 / coeffs[j]
				if c != c || math.IsInf(c, 0) {
					c = 1
				}
			}
			return c
		}

		if !causal && j > 0 && h > 0 {
			// TODO Compensate this.
			lookahead := int(math.Ceil(float64(n.hop)*c())) * 5
			future := j + lookahead
			filt := future-firstone > n.root.fs*n.root.opts.TransientMs/1000
			if future < len(lout) && coeffs[future] == 1 && filt {
				// Don't do anti-causal part if there is an overlap with the reset grain.
				if lookahead < n.nbuf {
					goto next
				}
				pin = j
				h = -n.hop
				j = future
				println(`future:`, pin, `-`, future, `, lah:`, lookahead)
			}
		}
	next:

		_ = firstone
		d := j - lastone
		if c() == 1 {
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
			copy(ingrain[ch][max(0, -i+n.nbuf/2):], input[ch][max(0, i-n.nbuf/2):clamp(0, len(lin), i+n.nbuf/2)])
		}

		if n.root.opts.Quality >= 0 && c() == 1 {
			if h > 0 {
				println(`reset past:`, j)
				n.resetPast(ingrain)
			} else if h < 0 {
				println(`reset future:`, j)
				n.resetFuture(ingrain)
			}
			n.bypassGrain(ingrain, grainbuf)
		} else {
			d := int(h / n.hop)
			if j == pin {
				d = 0
				println(`collision:`, pin)
			}
			n.advance([][][]float64{ingrain}, grainbuf, abs(c()), d)
		}

		// Cut pre-echo in transient regions.
		// TODO Works strange after implementing non-causality.
		if h > 0 && c() != 1 && d < n.nbuf/2 {
			for ch := range grainbuf {
				rr := grainbuf[ch][max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				// println(n.nbuf/2, d, n.hop, j, lastone)
				fill(grainbuf[ch][:n.nbuf/2-d-n.hop], 0)
			}
		}

		if n.root.opts.Onsets && c() != 1 {
			for ch := range grainbuf {
				clear(grainbuf[ch])
			}
		}

		// if h > 0 {
		loutgrain := lout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		add(loutgrain, grainbuf[0][clamp(0, n.nbuf, -j):])
		// }
		// if h < 0 {
		routgrain := rout[max(0, j-n.nbuf/2):clamp(0, len(lout), j+n.nbuf/2)]
		add(routgrain, grainbuf[1][clamp(0, n.nbuf, -j):])
		// }

		if h < 0 && j <= pin {
			h = n.hop
			j = lastone
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

func (n *warper) resetFuture(present [][]float64) {
	a := &n.a
	clear(a.Mid)
	for ch := range present {
		floats.Add(a.Mid, present[ch])
		n.enfft(a.C[ch], a.W, present[ch])
	}
	n.enfft(a.X, a.W, a.Mid)
	for w := range a.Ph {
		a.F[w], a.Future[w] = cmplx.Polar(a.X[w])
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
		// scale(output[ch], (1+n.wgain*float64(n.nbuf)/float64(n.hop))/2)
		// Where this is must be the last line:
		// scale(output[ch], n.wgain*float64(n.nbuf)/float64(n.hop))
		// Probably there is something wrong in gonum FFT implementation.
	}
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
// If direction is 1, it is doing forward integration.
// If it is -1, it is doing backward integration.
// If it is 0, integration is performed in both directions.
func (n *warper) advance(fs [][][]float64, output [][]float64, stretch float64, direction int) {
	a := &n.a

	present := fs[0]

	clear(a.Mid)
	for ch := range present {
		floats.Add(a.Mid, present[ch])
		n.enfft(a.C[ch], a.W, present[ch])
	}

	n.enfft(a.X, a.W, a.Mid)
	n.enfft(a.Xd, a.Wd, a.Mid)
	n.enfft(a.Xt, a.Wt, a.Mid)

	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	// TODO Works ONLY with 2x FFT oversampling.
	tadv := gettadv(a.X, a.Xd, float64(n.olap)*n.osamp)
	fadv := getfadv(a.X, a.Xt, stretch*2./n.osamp)

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
		a.E[w] = p
	}
	for ch := range present {
		for w := range a.X {
			a.C[ch][w] /= a.E[w]
		}
	}

	if direction == 0 {
		n.heap = n.heap[:n.nbins*2]
	} else {
		n.heap = n.heap[:n.nbins]
	}
	clear(n.heap)
	clear(n.arm)

	for j := range a.X {
		n.arm[j] = true
		if direction == 0 && a.P[j] > 1e-6 || direction > 0 {
			n.heap[j] = heaptriple{a.P[j], j, -1}
		}
		if direction == 0 && a.F[j] > 1e-6 || direction < 0 {
			n.heap[j] = heaptriple{a.F[j], j, 1}
		}
	}
	heapInit(&n.heap)

	// Do PGHI.
	for len(n.heap) > 0 {
		h := heapPop(&n.heap)
		w := h.w
		t := h.t
		// TODO Multiframe
		// https://ltfat.org/notes/ltfatnote040.pdf, page 5
		switch t {
		case 1:
			if n.arm[w] {
				adv := tadv(w)
				a.Ph[w] = a.Future[w] - adv
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, t - 1})
			}
		case -1:
			if n.arm[w] {
				adv := tadv(w)
				a.Ph[w] = a.Past[w] + adv
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, t + 1})
			}
		case 0:
			if w > 1 && n.arm[w-1] {
				adv := -a.Fadv[w-1]
				a.Ph[w-1] = a.Ph[w] + adv
				n.arm[w-1] = false
				heapPush(&n.heap, heaptriple{a.M[w-1], w - 1, t})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				adv := a.Fadv[w+1]
				a.Ph[w+1] = a.Ph[w] + adv
				n.arm[w+1] = false
				heapPush(&n.heap, heaptriple{a.M[w+1], w + 1, t})
			}
		}
	}

	if direction > 0 {
		copy(a.P, a.M)
	} else if direction < 0 {
		copy(a.F, a.M)
	}
	for w := range a.Ph {
		// Add stereo phase differences back through multiplication.
		a.E[w] = cmplx.Rect(1, a.Ph[w])
	}
	for w := range a.Ph {
		if direction > 0 {
			a.Past[w] = princarg(a.Ph[w])
		} else if direction < 0 {
			a.Future[w] = princarg(a.Ph[w])
		}
	}
	for ch := range output {
		cmplxs.MulTo(a.Co[ch], a.C[ch], a.E)
		n.defft(output[ch], a.Co[ch], true)
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
