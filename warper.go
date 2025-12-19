package nanowarp

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"

	"gonum.org/v1/gonum/dsp/fourier"
)

type warper struct {
	nfft    int
	nbuf    int
	nbins   int
	hop     int
	masking bool

	fft  *fourier.FFT
	arm  []bool
	norm float64
	heap hp
	img  [][]float64

	a wbufs
}
type wbufs struct {
	S, Mid, M, P, Phase, Pphase []float64 // Scratch buffers
	Fadv, Pfadv                 []float64
	Xprevs                      []complex128 // Lookahead framebuffer
	W, Wd, Wt                   []float64    // Window functions
	X, Xd, Xt, L, R, Lo, Ro     []complex128 // Complex spectra
}

func warperNew(nbuf int) (n *warper) {
	nextpow2 := int(math.Floor(math.Pow(2, math.Ceil(math.Log2(float64(nbuf))))))
	nfft := nextpow2 * 2
	nbins := nfft/2 + 1
	olap := 4
	n = &warper{
		nfft:  nfft,
		nbins: nbins,
		nbuf:  nbuf,
		hop:   nbuf / olap,
	}
	a := &n.a

	makeslices(a, nbins, nfft)

	// Exceptions.
	a.Phase = make([]float64, nbins)
	n.arm = make([]bool, nbins)

	hann(a.W[:nbuf])
	hannDx(a.Wd[:nbuf])
	windowT(a.W[:nbuf], a.Wt[:nbuf])
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
		} else {
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

	enfft(a.X, a.W, a.Mid)
	enfft(a.Xd, a.Wd, a.Mid)
	enfft(a.Xt, a.Wt, a.Mid)

	for w := range a.X {
		// Encode stereo phase differences and stretch mid only, keep original magnitudes.
		// NB: Phase difference in polar coordinates is complex division in cartesian.
		//     Phase sum is conversely a multiply.
		//     Hypot and multiplication are always cheaper than Atan2 and Sincos.
		//
		// See “Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a
		// real-time realization for a commercial music product.”
		m := mag(a.X[w])
		n := a.X[w] / complex(m, 0)
		if m < 1e-6 {
			n = complex(1, 0)
		}
		a.L[w] = a.L[w] / n
		a.R[w] = a.R[w] / n
		a.M[w] = m
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

	// Here we are using time-frequency reassignment[¹] as a way of obtaining
	// phase derivatives. Probably in future these derivatives will be replaced
	// with differences because currently phase is leaking in time domain and
	// finite differences guarrantee idempotency at stretch=1.
	//
	// [¹]: Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
	olap := float64(n.nbuf / n.hop)
	tadv := func(j int) float64 {
		if mag(a.X[j]) < 1e-6 {
			return 0
		}
		// TODO osampc is wrong. Does not work with non-power of 2 FFT sizes?
		osampc := float64(n.nfft/n.nbuf) / 2
		return (math.Pi*float64(j) + imag(a.Xd[j]/a.X[j])) / (olap * osampc)

	}
	fadv := getfadv(a.X[:n.nbins], a.Xt, stretch)
	for w := range a.X {
		a.Fadv[w] = princarg(fadv(w))
	}

	// Match phase reset points in stretched and original.
	for w := 1; w < n.nbins-1; w++ {
		m := func(j int) bool { return a.Pfadv[j]-a.Fadv[j] > math.Pi }
		if m(w) && m(w-1) && m(w+1) {
			a.Phase[w-1] = cmplx.Phase(a.X[w-1])
			a.Phase[w+1] = cmplx.Phase(a.X[w+1])
			a.Phase[w] = cmplx.Phase(a.X[w])
			n.arm[w] = false
			n.arm[w+1] = false
			n.arm[w-1] = false
		}
	}
	copy(a.Pfadv, a.Fadv)

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
		// Add stereo phase differences back through multiplication.
		a.Lo[w] = a.L[w] * cmplx.Rect(1, a.Phase[w])
		a.Ro[w] = a.R[w] * cmplx.Rect(1, a.Phase[w])
		// a.Ro[w] = cmplx.Rect(mag(a.R[w]), a.Phase[w])
		a.Pphase[w] = princarg(a.Phase[w])
	}

	defft := func(x []complex128, grain []float64) {
		n.fft.Sequence(a.S, x)
		for j := range a.S {
			a.S[j] /= n.norm
		}
		mul(a.S, a.W)
		copy(grain, a.S)
	}
	defft(a.Lo, loutgrain)
	defft(a.Ro, routgrain)
}
