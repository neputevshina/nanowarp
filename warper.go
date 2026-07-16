package nanowarp

import (
	"io"
	"math"
	"math/bits"
	"math/cmplx"
	"slices"

	"github.com/neputevshina/nanowarp/dspio"
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

	fft *fourier.FFT

	// PGHI-related
	arm      []bool    // Done mask
	heap     hp        // Heap
	arrows   [][2]int  // Integration directions for next frame
	parrows  [][2]int  // Current integration directions
	trace    []float64 // Ridge accumulator
	ftrace   []float64 // Filtered trace
	ridges   []uint    // Extracted ridges
	resetnow []bool    // Forced per-bin resets

	// PGPI-related
	reverse       [][]int      // Reverse index
	pairs, fpairs []heaptriple // Sorted tapes

	norm, wgain float64 // Global normalization factor and grain-only normalization factor

	a wbufs
}

type wbufs struct {
	S, Mid        []float64      // Scratch buffers
	F, M, P       []float64      `size:"nbins"` // Magnitudes: future, current and previous
	Ph            []float64      `size:"nbins"` // Current phase
	Past          []float64      // Phase accumulator
	Fadv, Tadv    []float64      // Partial derivatives
	W, Wr, Wd, Wt []float64      // Window functions
	X, Y, Xd, Xt  []complex128   // Complex spectra
	C, Co         [][]complex128 // Channel differences, original channels
}

func warperNew(nbuf, osamp, nch int, nanowarp *Nanowarp) (n *warper) {
	// FIXME Only 2x oversampling works, no more, no less.
	olap := 4 // Depends on window, see comment on blackmanHarrisClassic function.
	nfft := nextpow2(nbuf * osamp)
	n = &warper{
		nfft:  nfft,
		nbins: nfft/2 + 1,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		olap:  olap,
		osamp: float64(osamp),
		root:  nanowarp,
		lah:   3 * olap,
	}
	a := &n.a

	makeslices(a, n.nbins, nfft, nch, n.lah)
	n.arm = make([]bool, n.nbins)

	s := func(w []float64) []float64 {
		return w[:nbuf]
	}
	blackmanHarrisClassic(s(a.W))
	windowDx(s(a.Wd), s(a.W))
	windowT(s(a.Wt), s(a.W))
	n.wgain = windowDualUniform(s(a.Wr), s(a.W), n.hop)

	n.norm = float64(nfft) * float64(n.olap) * n.wgain

	n.heap = make(hp, n.nbins)
	n.fft = fourier.NewFFT(nfft)
	n.heap = make(hp, 2*n.nbins) // 2 for future and past.
	n.parrows = make([][2]int, 0, n.nbins)
	n.arrows = slices.Clone(n.parrows)
	n.ridges = make([]uint, n.nbins)
	n.trace = make([]float64, n.nbins)
	n.ftrace = make([]float64, n.nbins)

	return
}

func (n *warper) process6(in [][]float64, out [][]float64, phasor *Curve) {
	get := func() [][]float64 { return make2[float64](len(in), n.nfft) }
	nch := len(in)
	progress := n.root.opts.Progress

	lead := get()
	grain := get()
	crop := make([][]float64, nch)
	futurecrop := make([][]float64, nch)

	lastone := 0
	fivesec := n.root.fs * 5
	tsc := 0
	if progress != nil {
		progress <- Bp(0, 0)
	}
	for j := -n.nbuf / 2; j < len(out[0])-1+n.nbuf/2; j += n.hop {
		bounds := func(i int) int { return clamp(0, len(out[0])-1, i) }
		cut := func(s []float64, i int) []float64 { return s[bounds(i-n.nbuf/2):bounds(i+n.nbuf/2)] }
		i := int(phasor.ReverseSample(float64(j)))
		c := 1 / phasor.Dy(float64(j))

		if progress != nil && j/fivesec > tsc {
			progress <- Bp(float64(i), float64(j))
			tsc = j / fivesec
		}

		for ch := range nch {
			cr := cut(in[ch], i)
			if i < n.nbuf/2 {
				copy(lead[ch][max(0, n.nbuf/2-i):], cr)
				crop[ch] = lead[ch]
			} else {
				crop[ch] = cr
			}
		}

		q := n.root.opts.Resets
		normal, diff, _ := n.advance(crop, futurecrop, c, q >= -1 && c == 1, q == -1)
		n.synthesize(grain, normal, diff)

		d := j - lastone
		if c == 1 {
			lastone = j
		}
		for ch := range nch {
			// Cut pre-echo in transient regions.
			if c != 1 && d < n.nbuf/2 {
				rr := grain[ch][max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grain[ch][:n.nbuf/2-d-n.hop], 0)
			}

			if n.root.opts.Onsets && c != 1 {
				clear(grain[ch])
			}

			g := cut(out[ch], j)
			add(g, grain[ch][clamp(0, n.nbuf, -j):])
		}
	}
	if progress != nil {
		progress <- Bp(phasor.end.I, phasor.end.J)
		close(progress)
	}
}

func (n *warper) processFinal(in dspio.GrainSeeker, out dspio.GrainWriter, phasor *Curve) error {
	nch := in.NchRead()
	get := func() [][]float64 { return make2[float64](nch, n.nfft) }
	progress := n.root.opts.Progress

	lead := get()
	grain := get()
	crop := make([][]float64, nch)
	futurecrop := make([][]float64, nch)

	lastone := 0
	fivesec := n.root.fs * 5
	tsc := 0
	if progress != nil {
		progress <- Bp(0, 0)
	}
	var err error
	for j := -n.nbuf / 2; err == nil; j += n.hop {
		i := int(phasor.ReverseSample(float64(j)))
		c := 1 / phasor.Dy(float64(j))

		if progress != nil && j/fivesec > tsc {
			progress <- Bp(float64(i), float64(j))
			tsc = j / fivesec
		}

		err = in.GrainSeek(err, int64(i), lead)
		if err != nil && err != io.EOF {
			return err
		}

		q := n.root.opts.Resets
		normal, diff, _ := n.advance(crop, futurecrop, c, q >= -1 && c == 1, q == -1)
		n.synthesize(grain, normal, diff)

		d := j - lastone
		if c == 1 {
			lastone = j
		}
		for ch := range nch {
			// Cut pre-echo in transient regions.
			if c != 1 && d < n.nbuf/2 {
				rr := grain[ch][max(0, n.nbuf/2-d-n.hop) : n.nbuf/2-d]
				for i := range rr {
					rr[i] *= float64(i) / float64(len(rr))
				}
				fill(grain[ch][:n.nbuf/2-d-n.hop], 0)
			}

			if n.root.opts.Onsets && c != 1 {
				clear(grain[ch])
			}
		}
		_, err = out.SignalWrite(nil, grain)
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return err
		}
	}
	if progress != nil {
		progress <- Bp(phasor.end.I, phasor.end.J)
		close(progress)
	}
	return err
}

// advance constructs the next frame of the output.
func (n *warper) advance(ingrain, futuregrain [][]float64, stretch float64, reset, allreset bool) (normal []complex128, diff [][]complex128, mag []float64) {
	a := &n.a
	hp := n.root.opts.Hyperparams
	nch := len(ingrain)
	enfft := func(x []complex128, w, grain []float64) {
		clear(a.S)
		copy(a.S, grain)
		mul(a.S, w)
		n.fft.Coefficients(x, a.S)
	}

	clear(a.Mid)
	for ch := range ingrain {
		add(a.Mid, ingrain[ch])
		enfft(a.C[ch], a.W, ingrain[ch])
		copy(a.Co[ch], a.C[ch])
	}
	enfft(a.X, a.W, a.Mid)

	cmplxs.Abs(a.M, a.X)
	var arr [][2]int
	if n.root.opts.Quality == -1 {
		arr = n.bruteforcearrows(a.P, a.M, n.parrows, n.ridges)
	} else {
		arr = n.pghiarrows(a.P, a.M, n.parrows, n.ridges)
	}

	// Encode stereo phase differences and stretch mid only, keep original magnitudes.
	// NB: Phase difference in polar coordinates is complex division in cartesian.
	//     Phase sum is conversely a multiply.
	//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
	//
	// See Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
	// real-time realization for a commercial music product
	for w := range a.X {
		a.Y[w] = safediv(a.X[w], complex(a.M[w], 0))
		for ch := range nch {
			a.C[ch][w] = safediv(a.C[ch][w], a.Y[w])
		}
	}

	// Calculate partial derivatives of the phase.
	//
	// See Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms
	//
	// TODO Probably it will be more numerically stable to limit the phase accuum to
	// 0..1 and scale back to -π..π range at the poltocar conversion.
	// fadv must return 0..1 accordingly, simply defer π multiplication until the end.
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)
	for w := range a.X {
		a.Fadv[w] = princarg(fadv(a.X[:n.nbins], a.Xt, stretch, w))
		a.Tadv[w] = tadv(a.X[:n.nbins], a.Xd, float64(n.nfft/n.hop), w)
	}

	n.pghiintegrate(arr, a.Fadv, a.Tadv, a.Ph, a.Past)

	trace := n.trackridges(n.ftrace, n.trace, n.ridges, hp.HighRidgeHeight, hp.InfluenceRadius)

	// Bypass short ridges on phase reset.
	c := float64(hp.LongRidgeLength) * stretch
	upto := hztobin(hp.ResetUpToHz, n.nfft, n.root.fs)
	for w := range a.Y {
		if !reset || !allreset && trace[w] > c && w < upto {
			// Receive normals from the current phase, if not resetting.
			a.Y[w] = cmplx.Rect(1, a.Ph[w])
			continue
		}
		a.Ph[w] = cmplx.Phase(a.X[w])
		a.Y[w] = 1
		for ch := range nch {
			a.C[ch][w] = a.Co[ch][w]
		}
	}

	copy(a.P, a.M)
	copy(a.Past, a.Ph)

	return a.Y, a.C, a.M
}

func (n *warper) synthesize(outgrain [][]float64, normal []complex128, diff [][]complex128) {
	a := &n.a
	defft := func(out []float64, x []complex128) {
		n.fft.Sequence(a.S, x)
		floats.Scale(1/n.norm, a.S)
		mul(a.S, a.Wr)
		copy(out, a.S)
	}

	for ch := range outgrain {
		// Add stereo phase differences back through multiplication.
		mul(diff[ch], normal)
		defft(outgrain[ch], diff[ch])
	}
}

// Bitmasks for integration directions.
// Each ridge is encoded as a 3-bit mask.
const (
	right = 1 << iota
	down
	up
	ridgemask = right | down | up
)

func (n *warper) trackridges(out, trace []float64, ridges []uint, HighRidgeHeight, InfluenceRadius int) []float64 {
	for w, v := range ridges {
		p := boolfloat(bits.OnesCount(v&(ridgemask<<3)) >= 2)
		trace[w] = trace[w]*p + p
	}
	// Propagate vertically.
	l := -1
	for i := range trace {
		if l < 0 && trace[i] != 0 {
			l = i
		}
		if l >= 0 && trace[i] == 0 {
			v := slices.Max(trace[l:i])
			// Reset the track on a PGHI-detected transient.
			if i-l >= HighRidgeHeight {
				v = 0
			}
			fill(trace[l:i], v)
		}
		if trace[i] == 0 {
			l = -1
		}
	}
	if l > 0 {
		fill(trace[l:], slices.Max(trace[l:]))
	}
	// Propagate each trace to its native (per PGHI directions) region of influence,
	// limited by InfluenceRadius hyperparameter.
	clear(out)
	for w, v := range trace {
		if v == 0 {
			continue
		}
		for e := w - 1; e >= 0 && w-e <= InfluenceRadius; e-- {
			if trace[e] == 0 && ridges[e]&(down<<3) > 0 {
				out[e] = trace[w]
			} else {
				break
			}
		}
		for e := w + 1; e < len(trace) && e-w <= InfluenceRadius; e++ {
			if trace[e] == 0 && ridges[e]&(up<<3) > 0 {
				out[e] = trace[w]
			} else {
				break
			}
		}
	}
	// Add original traces to the output.
	for w, v := range trace {
		if v == 0 {
			continue
		}
		out[w] = v
	}
	return out
}

// pghiarrows calculates integration directions from a pair of magnitude
// buffers using priority queue.
//
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) pghiarrows(M, F []float64, arrows [][2]int, ridges []uint) [][2]int {
	n.heap = n.heap[:n.nbins]
	fill(n.arm, true)
	for w := range M {
		n.heap[w] = heaptriple{M[w], w, -1, 0}
		// Note that trajectories are calculated with a one frame lag.
		// Still good for our purposes.
		ridges[w] <<= 3
	}
	heapInit(&n.heap)
	arrows = arrows[:0]

	narm := n.nbins
	for len(n.heap) > 0 {
		if narm == 0 {
			break
		}
		h := heapPop(&n.heap)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[w] {
				arrows = append(arrows, [2]int{w, right})
				ridges[w] |= right << 3
				n.arm[w] = false
				narm--
				heapPush(&n.heap, heaptriple{F[w], w, 0, 0})
			}
		case 0:
			if w >= 1 && n.arm[w-1] {
				arrows = append(arrows, [2]int{w, down})
				ridges[w] |= down
				n.arm[w-1] = false
				narm--
				heapPush(&n.heap, heaptriple{F[w-1], w - 1, 0, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				arrows = append(arrows, [2]int{w, up})
				ridges[w] |= up
				n.arm[w+1] = false
				narm--
				heapPush(&n.heap, heaptriple{F[w+1], w + 1, 0, 0})
			}
		}
	}

	return arrows
}

// bruteforcearrows calculates integration directions from local magnitude maxima.
//
// It is faster than PGHI, but less accurate.
func (n *warper) bruteforcearrows(M, F []float64, arrows [][2]int, ridges []uint) [][2]int {
	fill(n.arm, true)
	for w := range M {
		// Note that trajectories are calculated with a one frame lag.
		// Still good for our purposes.
		ridges[w] <<= 3
	}
	arrows = arrows[:0]

	// Draw arrows from top to bottom, adding each from direction of largest neighbor.
	for w := n.nbins - 1; w >= 0; w-- {
		d := 0
		top := 0.
		if M[w] >= top {
			d = right
			top = M[w]
		}
		if w > 0 && F[w-1] > top {
			d = up
			top = F[w-1]
		}
		if w < n.nbins-1 && F[w+1] > top && n.arm[w+1] {
			d = down
		}

		switch d {
		case right:
			arrows = append(arrows, [2]int{w, right})
			ridges[w] |= right << 3
		case down:
			arrows = append(arrows, [2]int{w + 1, down})
			ridges[w+1] |= down
		case up:
			arrows = append(arrows, [2]int{w - 1, up})
			ridges[w-1] |= up
			// Make sure we don't add arrows pointing back.
			n.arm[w] = false
		}
	}

	// t := make([]float64, n.nbins)
	// for w := range arrows {
	// 	t[w] = float64(bits.OnesCount(n.ridges[w] & 0b111000))
	// }
	// oscope.Enable = true
	// oscope.Oscope(t, oscope.Name(`asdf`))

	// Repair ordering: start from rights, then do downs, then do ups reversed.
	do := func(arrows [][2]int, what, save int) int {
		arrows = arrows[save:]
		save = 0
		for i, a := range arrows {
			if a[1] == what {
				arrows[save], arrows[i] = arrows[i], arrows[save]
				save++
			}
		}
		return save
	}

	save := 0
	save += do(arrows, right, save)
	save += do(arrows, down, save)
	do(arrows, up, save)
	slices.Reverse(arrows[save:])

	return arrows
}

// pghiintegrate integrates partial derivatives of phase using directions
// obtained from pghiarrows.
func (n *warper) pghiintegrate(arrows [][2]int, Fadv, Tadv, Ph, Past []float64) {
	for _, e := range arrows {
		switch e[1] {
		case right:
			Ph[e[0]] = princarg(Past[e[0]]) + Tadv[e[0]]
		case down:
			Ph[e[0]-1] = Ph[e[0]] - Fadv[e[0]-1]
		case up:
			Ph[e[0]+1] = Ph[e[0]] + Fadv[e[0]+1]
		}
	}
}

// fadv calculates the partial derivative of the phase with respect
// to frequency using time-frequency reassignment.
func fadv(x, xt []complex128, stretch float64, w int) float64 {
	if cmplx.Abs(x[w]) == 0 {
		return 0
	}
	// NOTE Try len(x)-1 instead. Sounds worse on my $4 speakers.
	return -real(xt[w]/x[w])/float64(len(x))*math.Pi*stretch - math.Pi/2
}

// tadv calculates the partial derivative of the phase with respect
// to time using time-frequency reassignment.
func tadv(x, xd []complex128, olap float64, w int) float64 {
	if cmplx.Abs(x[w]) < 1e-6 {
		return 0
	}
	return (math.Pi*float64(w) + imag(xd[w]/x[w])) / (olap / 2)
}
