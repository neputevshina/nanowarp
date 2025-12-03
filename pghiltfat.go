package nanowarp

import (
	"math"
	"math/cmplx"
	"math/rand/v2"
)

type ltfat_int = int

func (n *Nanowarp) pghiltfat(stretch float64, out []complex128) {

}

type RTPGHIState struct {
	P       *RTPGHIUpdatePlan
	M, A, W int
	S       []float64
	TGrad   []float64
	FGrad   []float64
	Phase   []float64
	PhaseIn []float64
	Stretch float64
}

type RTPGHIUpdatePlan struct {
	H            []heaptriple
	DoneMask     []int
	Tol          float64
	M            int
	RandPhase    []float64
	RandPhaseLen int
	RandPhaseId  int
}

// AUDIT OFFBYONE
func shiftColsLeft(cols []float64, height, N int, newcol []float64) {
	for n := 0; n < N-1; n++ {
		copy(cols[n*height:(n+1)*height], cols[(n+1)*height:(n+2)*height])
	}
	if newcol != nil {
		copy(cols[(N-1)*height:], newcol[:height])
	} else {
		for i := (N - 1) * height; i < N*height; i++ {
			cols[i] = 0
		}
	}
}

func rtpghiAbs(in []complex128, height int, out []float64) {
	for i := 0; i < height; i++ {
		out[i] = cmplx.Abs(in[i])
	}
}

func rtpghiPhase(in []complex128, height int, out []float64) {
	for i := 0; i < height; i++ {
		out[i] = cmplx.Phase(in[i])
	}
}

func rtpghiGetStretch(p *RTPGHIState) float64 {
	return p.Stretch
}

func rtpghiInitStretch(p *RTPGHIState, stretch float64) {
	p.Stretch = stretch
}

func rtpghiSetTol(p *RTPGHIState, tol float64) error {
	if tol <= 0 || tol >= 1 {
		panic("tol must be in range ]0,1[")
	}
	p.P.Tol = tol
	return nil
}

func rtpghiInit(W, a, M int, tol float64) (*RTPGHIState, error) {
	if W <= 0 || a <= 0 || M <= 0 {
		panic("W, a, M must be positive")
	}
	if tol <= 0 || tol >= 1 {
		panic("tol must be in range ]0,1[")
	}

	M2 := M/2 + 1
	p := &RTPGHIState{}
	var err error
	p.P, err = rtpghiupdateInit(M, W, tol)
	if err != nil {
		return nil, err
	}
	p.S = make([]float64, 3*M2*W)
	p.TGrad = make([]float64, 2*M2*W)
	p.FGrad = make([]float64, 1*M2*W)
	p.Phase = make([]float64, 1*M2*W)
	p.PhaseIn = make([]float64, 3*M2*W)
	p.M = M
	p.A = a
	p.W = W
	p.Stretch = 1.0
	return p, nil
}

// AUDIT SUS
func rtpghiReset(p *RTPGHIState, sinit [][]float64) {
	M2 := p.M/2 + 1
	W := p.W
	clear(p.S)
	clear(p.TGrad)
	clear(p.FGrad)
	clear(p.Phase)
	clear(p.PhaseIn)
	if sinit != nil {
		for w := 0; w < W; w++ {
			if sinit[w] != nil {
				copy(p.S[2*w*M2:], sinit[w][:2*M2])
			}
		}
	}
}

func rtpghiExecute(p *RTPGHIState, cin []complex128, stretch float64, cout []complex128) {
	M2 := p.M/2 + 1
	W := p.W
	asyn := p.A
	aanaprev := int(math.Round(float64(asyn) / p.Stretch))
	aananext := int(math.Round(float64(asyn) / stretch))

	for w := 0; w < W; w++ {
		sCol := p.S[3*w*M2:]
		tgradCol := p.TGrad[2*w*M2:]
		fgradCol := p.FGrad[w*M2:]
		phaseCol := p.Phase[w*M2:]
		phaseinCol := p.PhaseIn[3*w*M2:]

		shiftColsLeft(sCol, M2, 3, nil)
		shiftColsLeft(tgradCol, M2, 2, nil)
		shiftColsLeft(phaseinCol, M2, 3, nil)

		rtpghiAbs(cin[w*M2:], M2, sCol[2*M2:])
		rtpghiPhase(cin[w*M2:], M2, phaseinCol[2*M2:])

		rtpghitgrad(phaseinCol, aanaprev, aananext, p.M, p.Stretch, tgradCol[M2:])

		if math.Abs(stretch-1.0) < 1e-4 {
			copy(phaseCol, phaseinCol[M2:M2+M2])
		} else {
			rtpghifgrad(phaseinCol[M2:], p.M, p.Stretch, fgradCol)
			rtpghiupdateExecute(p.P, sCol, tgradCol, fgradCol, phaseCol, phaseCol)
		}
		rtpghimagphase(sCol[M2:], phaseCol, M2, cout[w*M2:])
	}
	p.Stretch = stretch
}

func rtpghiDone(p **RTPGHIState) {
	if *p == nil {
		return
	}
	rtpghiupdateDone(&(*p).P)
	*p = nil
}

func rtpghiupdateInit(M, W int, tol float64) (*RTPGHIUpdatePlan, error) {
	M2 := M/2 + 1
	p := &RTPGHIUpdatePlan{}
	// AUDIT(aye): DoneMask is supposed to be two frames thick, no?
	p.DoneMask = make([]int, M2)
	p.RandPhaseLen = 10 * M2 * W
	p.RandPhase = make([]float64, p.RandPhaseLen)
	for i := range p.RandPhase {
		p.RandPhase[i] = 2.0 * math.Pi * rand.Float64()
	}
	p.Tol = tol
	p.M = M
	p.RandPhaseId = 0
	p.H = make([]heaptriple, 0, 2*M2)
	return p, nil
}

func rtpghiupdateExecuteWithMask(p *RTPGHIUpdatePlan, s, tgrad, fgrad, startphase []float64, mask []int, phase []float64) {
	copy(p.DoneMask, mask)
	rtpghiupdateExecuteCommon(p, s, tgrad, fgrad, startphase, phase)
}

func rtpghiupdateExecute(p *RTPGHIUpdatePlan, s, tgrad, fgrad, startphase, phase []float64) {
	// M2 := p.M/2 + 1
	for i := range p.DoneMask {
		p.DoneMask[i] = 0
	}
	rtpghiupdateExecuteCommon(p, s, tgrad, fgrad, startphase, phase)
}

// AUDIT(aye): Mostly LGTM
func rtpghiupdateExecuteCommon(p *RTPGHIUpdatePlan, s, tgrad, fgrad, startphase, phase []float64) {
	// h := p.H
	M2 := p.M/2 + 1
	quickbreak := M2
	oneover2 := 0.5
	donemask := p.DoneMask
	slog2 := s[M2:]

	logabstol := s[0]
	for m := 1; m < 2*M2; m++ {
		if s[m] > logabstol {
			logabstol = s[m]
		}
	}
	logabstol *= p.Tol

	p.H = p.H[:0] // Heap reset

	for m := 0; m < M2; m++ {
		if donemask[m] > 0 {
			heappush(&p.H, heaptriple{s[(m + M2)], (m + M2), 0})
			// heappush(&p.H, (m + M2))
			quickbreak--
		} else {
			if slog2[m] <= logabstol {
				donemask[m] = -1
				quickbreak--
			} else {
				heappush(&p.H, heaptriple{s[(m)], (m), 0})
				// heappush(&p.H, (m))
			}
		}
	}

	for quickbreak > 0 {
		// w := h.Delete()
		t := heappop(&p.H)
		w := t.w
		if w < 0 {
			break
		}
		if w >= M2 {
			wprev := w - M2
			if wprev != M2-1 && donemask[wprev+1] == 0 {
				phase[wprev+1] = phase[wprev] + (fgrad[wprev]+fgrad[wprev+1])*oneover2
				donemask[wprev+1] = 1
				// heappush(&p.H, (w + 1))
				heappush(&p.H, heaptriple{s[(w + 1)], (w + 1), 0})
				quickbreak--
			}
			if wprev != 0 && donemask[wprev-1] == 0 {
				phase[wprev-1] = phase[wprev] - (fgrad[wprev]+fgrad[wprev-1])*oneover2
				donemask[wprev-1] = 1
				// heappush(&p.H, (w - 1))
				heappush(&p.H, heaptriple{s[(w - 1)], (w - 1), 0})
				quickbreak--
			}
		} else {
			if donemask[w] == 0 {
				wnext := w + M2
				phase[w] = startphase[w] + (tgrad[w]+tgrad[wnext])*oneover2
				donemask[w] = 1
				heappush(&p.H, heaptriple{s[(wnext)], wnext, 0})
				// heappush(&p.H, (wnext))
				quickbreak--
			}
		}
	}

	for i := 0; i < M2; i++ {
		if donemask[i] < 0 {
			phase[i] = p.RandPhase[p.RandPhaseId]
			p.RandPhaseId++
			if p.RandPhaseId >= p.RandPhaseLen {
				p.RandPhaseId = 0
			}
		}
	}
}

func rtpghiupdateDone(p **RTPGHIUpdatePlan) {
	if *p == nil {
		return
	}
	(*p).H = nil
	*p = nil
}

func rtpghifgrad(phase []float64, M int, stretch float64, fgrad []float64) {
	M2 := M/2 + 1
	for m := 1; m < M2-1; m++ {
		fgrad[m] = (princarg(phase[m+1]-phase[m]) + princarg(phase[m]-phase[m-1])) / 2.0 * stretch
	}
	fgrad[0] = phase[0] * stretch
	fgrad[M2-1] = phase[M2-1] * stretch
}

func rtpghitgrad(phase []float64, aanaprev, aananext, M int, stretch float64, tgrad []float64) {
	M2 := M/2 + 1
	asyn := float64(aanaprev) * stretch
	const1prev := 2.0 * math.Pi * float64(aanaprev) / float64(M)
	const1next := 2.0 * math.Pi * float64(aananext) / float64(M)
	const2 := 2.0 * math.Pi * asyn / float64(M)
	pcol0 := phase
	pcol1 := phase[M2:]
	pcol2 := phase[2*M2:]

	for m := 0; m < M2; m++ {
		tgrad[m] = asyn*(princarg(pcol2[m]-pcol1[m]-const1next*float64(m))/(2.0*float64(aananext))-
			princarg(pcol1[m]-pcol0[m]-const1prev*float64(m))/(2.0*float64(aanaprev))) + const2*float64(m)
	}
}

func rtpghimagphase(s, phase []float64, L int, c []complex128) {
	for l := 0; l < L; l++ {
		c[l] = complex(s[l]*math.Cos(phase[l]), s[l]*math.Sin(phase[l]))
	}
}
