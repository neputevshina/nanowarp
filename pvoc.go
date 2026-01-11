package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"slices"

	"github.com/neputevshina/nanowarp/oscope"
	"gonum.org/v1/gonum/dsp/fourier"
)

type warper struct {
	nfft    int
	nbuf    int
	nbins   int
	hop     int
	masking bool
	diffadv bool
	root    *Nanowarp

	fft  *fourier.FFT
	arm  []bool
	norm float64
	heap hp
	img  [][]float64

	a wbufs
}
type wbufs struct {
	S, Mid, M, P, Phase, Pphase     []float64 // Scratch buffers
	Fadv, Pfadv, Ldiff, Rdiff       []float64
	Xprevs                          []complex128 // Lookahead framebuffer
	W, Wr, Wd, Wt                   []float64    // Window functions
	X, Xd, Xt, L, Ld, R, Rd, Lo, Ro []complex128 // Complex spectra
}

func warperNew(nbuf int, nanowarp *Nanowarp) (n *warper) {
	nextpow2 := int(math.Floor(math.Pow(2, math.Ceil(math.Log2(float64(nbuf))))))
	nfft := nextpow2 * 2
	nbins := nfft/2 + 1
	olap := 4
	n = &warper{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
		root:  nanowarp,
	}
	a := &n.a

	makeslices(a, nbins, nfft)

	// Exceptions.
	a.Phase = make([]float64, nbins)
	n.arm = make([]bool, nbins)

	niemitalo(a.W[:nbuf])
	windowDx(a.W[:nbuf], a.Wd[:nbuf])
	windowT(a.W[:nbuf], a.Wt[:nbuf])
	copy(a.Wr[:nbuf], a.W[:nbuf])
	slices.Reverse(a.Wr[:nbuf])
	n.norm = float64(nfft) / float64(n.hop) * float64(nfft) * windowGain(n.a.W)

	n.fft = fourier.NewFFT(nfft)

	return
}

func (n *warper) process(lin, rin, lout, rout []float64, stretch float64, delay float64) {
	inhop := float64(n.hop) / stretch
	ih, fh := math.Modf(inhop)
	fmt.Fprintln(os.Stderr, `(*warper).process`)
	fmt.Fprintln(os.Stderr, `stretch:`, stretch, `nbuf:`, n.nbuf, `nsampin:`, len(lin), `nsampout:`, len(lout))
	fmt.Fprintln(os.Stderr, `inhop:`, inhop, `whole:`, ih, `frac:`, fh, `interval:`, float64(n.hop)/(ih+1), `-`, float64(n.hop)/(ih))
	fmt.Fprintln(os.Stderr, `outhop:`, n.hop)

	n.start(lin, lout)
	n.start(rin, rout)
	lgrainbuf := make([]float64, n.nfft)
	rgrainbuf := make([]float64, n.nfft)

	id, fd := math.Modf(delay)
	j := n.hop + int(id)
	dh := fh
	dd := fd
	for i := int(ih); i < len(lin); i += int(ih) {
		if j > len(lout) {
			break
		}
		lingrain := lin[i:min(len(lin), i+n.nbuf)]
		ringrain := rin[i:min(len(lin), i+n.nbuf)]
		loutgrain := lout[j:min(len(lout), j+n.nbuf)]
		routgrain := rout[j:min(len(lout), j+n.nbuf)]

		// Dither both input and output hop sizes to get
		// fractional stretch.
		// Snap stretch coefficient to the dithered input hop size
		// to hopefully get more accurate phase.
		dh += fh
		hop := float64(n.hop) / (ih)
		if dh > 0 {
			i += 1
			dh -= 1
			hop = float64(n.hop) / (ih + 1)
		}
		n.advance(lingrain, ringrain, lgrainbuf, rgrainbuf, hop)
		dd += fd
		if dh > 0 {
			j += 1
			dd -= 1
		}
		add(loutgrain, lgrainbuf)
		add(routgrain, rgrainbuf)
		j += n.hop
	}
}

func (n *warper) process1(lin, rin, lout, rout []float64, onsets, stretch []float64, align int, stretchout []float64) {
	fmt.Fprintln(os.Stderr, `(*warper).process1`)
	fmt.Fprintln(os.Stderr, `outhop:`, n.hop)

	n.start(lin, lout)
	n.start(rin, rout)
	lgrainbuf := make([]float64, n.nfft)
	rgrainbuf := make([]float64, n.nfft)

	delay := func(stretch float64) float64 {
		return float64(align) * (stretch - 1) * 2
	}

	hop := float64(n.hop)
	ih, dh, fh := 0., 0., 0.
	id, fd := math.Modf(delay(stretch[0]))
	j := n.hop + int(id)
	dd := fd
	k := 0.

	for i := 0; i < len(lin); i += int(ih) {
		s := max(0, (stretch[i]*float64(i)-float64(j)+stretch[i]*4096)/4096)
		inhop := hop / s
		// Add 2 hop sizes to center the window.
		if onsets[clamp(0, len(onsets)-1, i+2*n.hop)] == 0 {
			inhop = hop
		}
		ih, fh = math.Modf(inhop)
		id, fd = math.Modf(delay(stretch[i]))

		if i > len(lin) {
			break
		}
		if j > len(lout) {
			break
		}
		lingrain := lin[i:min(len(lin), i+n.nbuf)]
		ringrain := rin[i:min(len(lin), i+n.nbuf)]
		loutgrain := lout[j:min(len(lout), j+n.nbuf)]
		routgrain := rout[j:min(len(lout), j+n.nbuf)]

		dh += fh
		h := hop / ih
		_ = h
		if dh > 0 {
			i += 1
			dh -= 1
			h = hop / (ih + 1)
		}
		n.advance(lingrain, ringrain, lgrainbuf, rgrainbuf, h)
		dd += fd
		if dh > 0 {
			j += 1
			dd -= 1
		}

		add(loutgrain, lgrainbuf)
		add(routgrain, rgrainbuf)

		k += hop / stretch[i]
		j += n.hop
	}
}

func (n *warper) start(in []float64, out []float64) {
	a := &n.a

	clear(a.S)
	copy(a.S, in[:min(len(in), n.nbuf)])
	mul(a.S, a.W)
	n.fft.Coefficients(a.X, a.S)

	for w := range a.X {
		a.M[w] = mag(a.X[w])
		// FIXME Right phase (not sum) is the starting previous phase.
		a.Pphase[w] = cmplx.Phase(a.X[w])
	}
	copy(a.P, a.M)

	for j := range a.S[:n.nbuf] {
		a.S[j] = in[:min(len(in), n.nbuf)][j] * a.W[j] * a.W[j] / n.norm * float64(n.nfft)
	}
	add(out[:min(len(out), n.nbuf)], a.S)
}

// advance adds to the phase of the output by one frame using
// phase gradient heap integration.
// See Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.
// (https://arxiv.org/pdf/2202.07382)
func (n *warper) advance(lingrain, ringrain, loutgrain, routgrain []float64, stretch float64) {
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
	if n.diffadv {
		enfft(a.Ld, a.Wd, lingrain)
		enfft(a.Rd, a.Wd, ringrain)
	}

	enfft(a.X, a.W, a.Mid)
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	// Here we are using time-frequency reassignment[¹] as a way of obtaining
	// phase derivatives. Probably in future these derivatives will be replaced
	// with differences because finite differences guarrantee idempotency at stretch=1.
	//
	// See also https://github.com/y-fujii/mini_pvdr/, which uses finite differences and
	// achieves less phase noise than current method. It both sounds clearer and on the
	// reassigned point density spectrogram there are almost no ”squiggly” trajectories.
	//
	// Elastiqué, however, behaves similarly to the current Nanowarp, but it is a 20 year old algorithm.
	//
	// [¹]: Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	olap := float64(n.nbuf / n.hop)
	osampc := float64(n.nfft / n.nbuf) // TODO osampc is wrong. Does not work with zero padded non-power of 2 FFT sizes?
	tadv := gettadv(a.X[:n.nbins], a.Xd, osampc, olap)
	ltadv := gettadv(a.L[:n.nbins], a.Ld, osampc, olap)
	rtadv := gettadv(a.R[:n.nbins], a.Rd, osampc, olap)
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch)

	// Force reset on a no-op stretch, including detected transients.
	if stretch == 1 {
		for w := range a.Phase {
			a.Pphase[w] = cmplx.Phase(a.X[w])
		}
		copy(a.P, a.M)
		copy(a.Lo, a.L)
		copy(a.Ro, a.R)
		goto skip
	}

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
		if n.diffadv {
			a.Ldiff[w] = tadv(w) - ltadv(w)
			a.Rdiff[w] = tadv(w) - rtadv(w)
		} else {
			a.L[w] /= p
			a.R[w] /= p
		}
	}

	n.heap = make(hp, n.nbins)
	clear(n.arm)

	for j := range a.X {
		n.arm[j] = true
		if n.masking {
			// Allow time-phase propagation only for local maxima.
			// Which is the simplest possible auditory masking model.
			if j == 0 || j == n.nbins-1 ||
				a.M[j-1] < a.M[j] && a.M[j+1] < a.M[j] {
				n.heap[j] = heaptriple{a.P[j], j, -1}
			}
		} else {
			n.heap[j] = heaptriple{a.P[j], j, -1}
		}
	}
	heapInit(&n.heap)

	if !n.root.opts.Noreset {
		// Match phase reset points in stretched and original.
		var dif []float64
		if oscope.Enable {
			dif = make([]float64, n.nbins)
		}
		vertical := 0
		for w := 1; w < n.nbins-1; w++ {
			m := func(j int) bool { return a.Pfadv[j]-a.Fadv[j] > math.Pi*3/2 }
			if m(w) && m(w-1) && m(w+1) {
				a.Phase[w-1] = cmplx.Phase(a.X[w-1])
				a.Phase[w+1] = cmplx.Phase(a.X[w+1])
				a.Phase[w] = cmplx.Phase(a.X[w])
				n.arm[w] = false
				n.arm[w+1] = false
				n.arm[w-1] = false
				if oscope.Enable {
					dif[w] = 1
					dif[w+1] = 1
					dif[w-1] = 1
				}
				vertical++
			}
		}
		copy(a.Pfadv, a.Fadv)
		// if vertical > n.nbins/4 {
		// 	for w := range a.Phase {
		// 		a.Pphase[w] = cmplx.Phase(a.X[w])
		// 	}
		// 	copy(a.P, a.M)
		// 	copy(a.Lo, a.L)
		// 	copy(a.Ro, a.R)
		// }
		oscope.Oscope(dif)
	}

	for len(n.heap) > 0 {
		h := heapPop(&n.heap).(heaptriple)
		w := h.w
		switch h.t {
		case -1:
			if n.arm[w] {
				a.Phase[w] = a.Pphase[w] + tadv(w)
				n.arm[w] = false
				heapPush(&n.heap, heaptriple{a.M[w], w, 0})
			}
		case 0:
			if w > 1 && n.arm[w-1] {
				a.Phase[w-1] = a.Phase[w] - a.Fadv[w-1]
				n.arm[w-1] = false
				heapPush(&n.heap, heaptriple{a.M[w-1], w - 1, 0})
			}
			if w < n.nbins-1 && n.arm[w+1] {
				a.Phase[w+1] = a.Phase[w] + a.Fadv[w+1]
				n.arm[w+1] = false
				heapPush(&n.heap, heaptriple{a.M[w+1], w + 1, 0})
			}
		}
	}

	copy(a.P, a.M)
	for w := range a.Phase {
		if n.diffadv {
			a.Lo[w] = cmplx.Rect(mag(a.L[w]), a.Phase[w]+a.Ldiff[w])
			a.Ro[w] = cmplx.Rect(mag(a.R[w]), a.Phase[w]+a.Rdiff[w])
		} else {
			// Add stereo phase differences back through multiplication.
			a.Lo[w] = a.L[w] * cmplx.Rect(1, a.Phase[w])
			a.Ro[w] = a.R[w] * cmplx.Rect(1, a.Phase[w])
		}
		a.Pphase[w] = princarg(a.Phase[w])
	}

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
