//go:build goexperiment.simd

package pffft

import (
	"simd"
)

// v4sf tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br);
// ar=VSUB(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VADD(ai,tmp); }
func vcplxmul(ar, ai, br, bi simd.Float32s) (simd.Float32s, simd.Float32s) {
	ar = ar.Mul(br).Sub(ai.Mul(bi))
	ai = ai.Mul(br).Add(ar.Mul(bi))
	return ar, ai
}

// v4sf tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br);
// ar=VADD(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VSUB(ai,tmp);
func vcplxmulconj(ar, ai, br, bi simd.Float32s) (simd.Float32s, simd.Float32s) {
	ar = ar.Mul(br).Add(ai.Mul(bi))
	ai = ai.Mul(br).Sub(ar.Mul(bi))
	return ar, ai
}

func svmul(s float32, ps simd.Float32s) simd.Float32s {
	return simd.BroadcastFloat32s(s).Mul(ps)
}

func vadd(a, b simd.Float32s) simd.Float32s {
	return a.Add(b)
}

func vsub(a, b simd.Float32s) simd.Float32s {
	return a.Sub(b)
}

func vmul(a, b simd.Float32s) simd.Float32s {
	return a.Mul(b)
}

func vmadd(a, b, c simd.Float32s) simd.Float32s {
	return a.MulAdd(b, c)
}

func ld_ps1(s float32) simd.Float32s {
	return simd.BroadcastFloat32s(s)
}

//go:noinline
func passf2_ps(ido, l1 int, cc, ch []simd.Float32s, wa1 []float32, fsign float32) {
	l1ido := l1 * ido
	if ido <= 2 {
		for k := 0; k < l1ido; k, ch, cc = k+ido, ch[ido:], cc[ido*2:] {
			ch[0] = vadd(cc[0], cc[ido+0])
			ch[l1ido] = vsub(cc[0], cc[ido+0])
			ch[1] = vadd(cc[1], cc[ido+1])
			ch[l1ido+1] = vsub(cc[1], cc[ido+1])
		}
	} else {
		for k := 0; k < l1ido; k, ch, cc = k+ido, ch[ido:], cc[ido*2:] {
			for i := 0; i < ido-1; i += 2 {
				tr2 := vsub(cc[i+0], cc[i+ido+0])
				ti2 := vsub(cc[i+1], cc[i+ido+1])
				wr := ld_ps1(wa1[i])
				wi := vmul(ld_ps1(fsign), ld_ps1(wa1[i+1]))
				ch[i] = vadd(cc[i+0], cc[i+ido+0])
				ch[i+1] = vadd(cc[i+1], cc[i+ido+1])
				tr2, ti2 = vcplxmul(tr2, ti2, wr, wi)
				ch[i+l1ido] = tr2
				ch[i+l1ido+1] = ti2
			}
		}
	}
}

//go:noinline
func passf3_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2 []float32, fsign float32) {
	const taur = -0.5
	taui := 0.866025403784439 * fsign
	if ido <= 2 {
		panic(`fatal: passf3_ps: ido <= 2`)
	}
	l1ido := l1 * ido
	for k := 0; k < l1ido; k, ch, cc = k+ido, ch[ido:], cc[ido*3:] {
		for i := 0; i < ido-1; i += 2 {
			tr2 := vadd(cc[i+ido], cc[i+2*ido])
			cr2 := vadd(cc[i], svmul(taur, tr2))
			ch[i] = vadd(cc[i], tr2)
			ti2 := vadd(cc[i+ido+1], cc[i+2*ido+1])
			ci2 := vadd(cc[i+1], svmul(taur, ti2))
			ch[i+1] = vadd(cc[i+1], ti2)
			cr3 := svmul(taui, vsub(cc[i+ido], cc[i+2*ido]))
			ci3 := svmul(taui, vsub(cc[i+ido+1], cc[i+2*ido+1]))
			dr2 := vsub(cr2, ci3)
			dr3 := vadd(cr2, ci3)
			di2 := vadd(ci2, cr3)
			di3 := vsub(ci2, cr3)
			wr1 := wa1[i]
			wi1 := fsign * wa1[i+1]
			wr2 := wa2[i]
			wi2 := fsign * wa2[i+1]
			dr2, di2 = vcplxmul(dr2, di2, ld_ps1(wr1), ld_ps1(wi1))
			ch[i+l1ido] = dr2
			ch[i+l1ido+1] = di2
			dr3, di3 = vcplxmul(dr3, di3, ld_ps1(wr2), ld_ps1(wi2))
			ch[i+2*l1ido] = dr3
			ch[i+2*l1ido+1] = di3

		}
	}
}

//go:noinline
func passf4_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2, wa3 []float32, fsign float32) {
	l1ido := l1 * ido
	if ido == 2 {
		for k := 0; k < l1ido; k, ch, cc = k+ido, ch[ido:], cc[ido*4:] {
			tr1 := vsub(cc[0], cc[2*ido+0])
			tr2 := vadd(cc[0], cc[2*ido+0])
			ti1 := vsub(cc[1], cc[2*ido+1])
			ti2 := vadd(cc[1], cc[2*ido+1])
			ti4 := vmul(vsub(cc[1*ido+0], cc[3*ido+0]), ld_ps1(fsign))
			tr4 := vmul(vsub(cc[3*ido+1], cc[1*ido+1]), ld_ps1(fsign))
			tr3 := vadd(cc[ido+0], cc[3*ido+0])
			ti3 := vadd(cc[ido+1], cc[3*ido+1])

			ch[0*l1ido+0] = vadd(tr2, tr3)
			ch[0*l1ido+1] = vadd(ti2, ti3)
			ch[1*l1ido+0] = vadd(tr1, tr4)
			ch[1*l1ido+1] = vadd(ti1, ti4)
			ch[2*l1ido+0] = vsub(tr2, tr3)
			ch[2*l1ido+1] = vsub(ti2, ti3)
			ch[3*l1ido+0] = vsub(tr1, tr4)
			ch[3*l1ido+1] = vsub(ti1, ti4)
		}
	} else {
		for k := 0; k < l1ido; k, ch, cc = k+ido, ch[ido:], cc[ido*4:] {
			for i := 0; i < ido-1; i += 2 {
				tr1 := vsub(cc[i+0], cc[i+2*ido+0])
				tr2 := vadd(cc[i+0], cc[i+2*ido+0])
				ti1 := vsub(cc[i+1], cc[i+2*ido+1])
				ti2 := vadd(cc[i+1], cc[i+2*ido+1])
				tr4 := vmul(vsub(cc[i+3*ido+1], cc[i+1*ido+1]), ld_ps1(fsign))
				ti4 := vmul(vsub(cc[i+1*ido+0], cc[i+3*ido+0]), ld_ps1(fsign))
				tr3 := vadd(cc[i+ido+0], cc[i+3*ido+0])
				ti3 := vadd(cc[i+ido+1], cc[i+3*ido+1])

				ch[i] = vadd(tr2, tr3)
				cr3 := vsub(tr2, tr3)
				ch[i+1] = vadd(ti2, ti3)
				ci3 := vsub(ti2, ti3)

				cr2 := vadd(tr1, tr4)
				cr4 := vsub(tr1, tr4)
				ci2 := vadd(ti1, ti4)
				ci4 := vsub(ti1, ti4)
				wr1 := wa1[i]
				wi1 := fsign * wa1[i+1]
				cr2, ci2 = vcplxmul(cr2, ci2, ld_ps1(wr1), ld_ps1(wi1))
				wr2 := wa2[i]
				wi2 := fsign * wa2[i+1]
				ch[i+l1ido] = cr2
				ch[i+l1ido+1] = ci2

				cr3, ci3 = vcplxmul(cr3, ci3, ld_ps1(wr2), ld_ps1(wi2))
				wr3 := wa3[i]
				wi3 := fsign * wa3[i+1]
				ch[i+2*l1ido] = cr3
				ch[i+2*l1ido+1] = ci3

				cr4, ci4 = vcplxmul(cr4, ci4, ld_ps1(wr3), ld_ps1(wi3))
				ch[i+3*l1ido] = cr4
				ch[i+3*l1ido+1] = ci4
			}
		}
	}
}

//go:noinline
func passf5_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2, wa3, wa4 []float32, fsign float32) {
	const tr11 = .309016994374947
	ti11 := .951056516295154 * fsign
	const tr12 = -.809016994374947
	ti12 := .587785252292473 * fsign

	if ido <= 2 {
		panic(`fatal: passf5_ps: ido <= 2`)
	}
	cc_ref := func(a_1, a_2 int) simd.Float32s { return cc[(a_2-1)*ido+a_1+1] }
	ch_ref := func(a_1, a_3 int) *simd.Float32s { return &ch[(a_3-1)*l1*ido+a_1+1] }

	for k := 0; k < l1; k, ch, cc = k+ido, ch[ido:], cc[ido*5:] {
		for i := 0; i < ido-1; i += 2 {
			ti5 := vsub(cc_ref(i, 2), cc_ref(i, 5))
			ti2 := vadd(cc_ref(i, 2), cc_ref(i, 5))
			ti4 := vsub(cc_ref(i, 3), cc_ref(i, 4))
			ti3 := vadd(cc_ref(i, 3), cc_ref(i, 4))
			tr5 := vsub(cc_ref(i-1, 2), cc_ref(i-1, 5))
			tr2 := vadd(cc_ref(i-1, 2), cc_ref(i-1, 5))
			tr4 := vsub(cc_ref(i-1, 3), cc_ref(i-1, 4))
			tr3 := vadd(cc_ref(i-1, 3), cc_ref(i-1, 4))
			*ch_ref(i-1, 1) = vadd(cc_ref(i-1, 1), vadd(tr2, tr3))
			*ch_ref(i, 1) = vadd(cc_ref(i, 1), vadd(ti2, ti3))
			cr2 := vadd(cc_ref(i-1, 1), vadd(svmul(tr11, tr2), svmul(tr12, tr3)))
			ci2 := vadd(cc_ref(i, 1), vadd(svmul(tr11, ti2), svmul(tr12, ti3)))
			cr3 := vadd(cc_ref(i-1, 1), vadd(svmul(tr12, tr2), svmul(tr11, tr3)))
			ci3 := vadd(cc_ref(i, 1), vadd(svmul(tr12, ti2), svmul(tr11, ti3)))
			cr5 := vadd(svmul(ti11, tr5), svmul(ti12, tr4))
			ci5 := vadd(svmul(ti11, ti5), svmul(ti12, ti4))
			cr4 := vsub(svmul(ti12, tr5), svmul(ti11, tr4))
			ci4 := vsub(svmul(ti12, ti5), svmul(ti11, ti4))
			dr3 := vsub(cr3, ci4)
			dr4 := vadd(cr3, ci4)
			di3 := vadd(ci3, cr4)
			di4 := vsub(ci3, cr4)
			dr5 := vadd(cr2, ci5)
			dr2 := vsub(cr2, ci5)
			di5 := vsub(ci2, cr5)
			di2 := vadd(ci2, cr5)
			wr1 := wa1[i]
			wi1 := fsign * wa1[i+1]
			wr2 := wa2[i]
			wi2 := fsign * wa2[i+1]
			wr3 := wa3[i]
			wi3 := fsign * wa3[i+1]
			wr4 := wa4[i]
			wi4 := fsign * wa4[i+1]
			dr2, di2 = vcplxmul(dr2, di2, ld_ps1(wr1), ld_ps1(wi1))
			*ch_ref(i-1, 2) = dr2
			*ch_ref(i, 2) = di2
			dr3, di3 = vcplxmul(dr3, di3, ld_ps1(wr2), ld_ps1(wi2))
			*ch_ref(i-1, 3) = dr3
			*ch_ref(i, 3) = di3
			dr4, di4 = vcplxmul(dr4, di4, ld_ps1(wr3), ld_ps1(wi3))
			*ch_ref(i-1, 4) = dr4
			*ch_ref(i, 4) = di4
			dr5, di5 = vcplxmul(dr5, di5, ld_ps1(wr4), ld_ps1(wi4))
			*ch_ref(i-1, 5) = dr5
			*ch_ref(i, 5) = di5
		}
	}
}

//go:noinline
func radf2_ps(ido, l1 int, cc, ch []simd.Float32s, wa1 []float32) {
	l1ido := l1 * ido
	for k := 0; k < l1ido; k += ido {
		a, b := cc[k], cc[k+l1ido]
		ch[2*k] = vadd(a, b)
		ch[2*(k+ido)-1] = vsub(a, b)
	}
	if ido < 2 {
		return
	}
	if ido != 2 {
		for k := 0; k < l1ido; k += ido {
			for i := 2; i < ido; i += 2 {
				tr2, ti2 := cc[i-1+k+l1ido], cc[i+k+l1ido]
				br, bi := cc[i-1+k], cc[i+k]
				tr2, ti2 = vcplxmulconj(tr2, ti2, ld_ps1(wa1[i-2]), ld_ps1(wa1[i-1]))
				ch[i+2*k] = vadd(bi, ti2)
				ch[2*(k+ido)-i] = vsub(ti2, bi)
				ch[i-1+2*k] = vadd(br, tr2)
				ch[2*(k+ido)-i-1] = vsub(br, tr2)
			}
		}
		if ido%2 == 1 {
			return
		}
	}
	for k := 0; k < l1ido; k += ido {
		ch[2*k+ido] = svmul(-1, cc[ido-1+k+l1ido])
		ch[2*k+ido-1] = cc[k+ido-1]
	}
}

//go:noinline
func radb2_ps(ido, l1 int, cc, ch []simd.Float32s, wa1 []float32) {
	l1ido := l1 * ido
	for k := 0; k < l1ido; k += ido {
		a, b := cc[2*k], cc[2*(k+ido)-1]
		ch[k] = vadd(a, b)
		ch[k+l1ido] = vsub(a, b)
	}
	if ido < 2 {
		return
	}
	if ido != 2 {
		for k := 0; k < l1ido; k += ido {
			for i := 2; i < ido; i += 2 {
				a, b := cc[i-1+2*k], cc[2*(k+ido)-i-1]
				c, d := cc[i+0+2*k], cc[2*(k+ido)-i+0]
				ch[i-1+k] = vadd(a, b)
				tr2 := vsub(a, b)
				ch[i+0+k] = vsub(c, d)
				ti2 := vadd(c, d)
				tr2, ti2 = vcplxmul(tr2, ti2, ld_ps1(wa1[i-2]), ld_ps1(wa1[i-1]))
				ch[i-1+k+l1ido] = tr2
				ch[i+0+k+l1ido] = ti2
			}
		}
		if ido%2 == 1 {
			return
		}
	}
	for k := 0; k < l1ido; k += ido {
		a := cc[2*k+ido-1]
		b := cc[2*k+ido]
		ch[k+ido-1] = vadd(a, a)
		ch[k+ido-1+l1ido] = svmul(-2, b)
	}
}

//go:noinline
func radf3_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2 []float32) {
	const taur = -0.5
	const taui = 0.866025403784439
	for k := range l1 {
		cr2 := vadd(cc[(k+l1)*ido], cc[(k+2*l1)*ido])
		ch[3*k*ido] = vadd(cc[k*ido], cr2)
		ch[(3*k+2)*ido] = svmul(taui, vsub(cc[(k+l1*2)*ido], cc[(k+l1)*ido]))
		ch[ido-1+(3*k+1)*ido] = vadd(cc[k*ido], svmul(taur, cr2))
	}
	if ido == 1 {
		return
	}
	for k := range l1 {
		for i := 2; i < ido; i += 2 {
			ic := ido - i
			wr1 := ld_ps1(wa1[i-2])
			wi1 := ld_ps1(wa1[i-1])
			dr2 := cc[i-1+(k+l1)*ido]
			di2 := cc[i+(k+l1)*ido]
			dr2, di2 = vcplxmulconj(dr2, di2, wr1, wi1)

			wr2 := ld_ps1(wa2[i-2])
			wi2 := ld_ps1(wa2[i-1])
			dr3 := cc[i-1+(k+l1*2)*ido]
			di3 := cc[i+(k+l1*2)*ido]
			dr3, di3 = vcplxmulconj(dr3, di3, wr2, wi2)

			cr2 := vadd(dr2, dr3)
			ci2 := vadd(di2, di3)
			ch[i-1+3*k*ido] = vadd(cc[i-1+k*ido], cr2)
			ch[i+3*k*ido] = vadd(cc[i+k*ido], ci2)
			tr2 := vadd(cc[i-1+k*ido], svmul(taur, cr2))
			ti2 := vadd(cc[i+k*ido], svmul(taur, ci2))
			tr3 := svmul(taui, vsub(di2, di3))
			ti3 := svmul(taui, vsub(dr3, dr2))
			ch[i-1+(3*k+2)*ido] = vadd(tr2, tr3)
			ch[ic-1+(3*k+1)*ido] = vsub(tr2, tr3)
			ch[i+(3*k+2)*ido] = vadd(ti2, ti3)
			ch[ic+(3*k+1)*ido] = vsub(ti3, ti2)
		}
	}
}

//go:noinline
func radb3_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2 []float32) {
	const taur = -0.5
	const taui = 0.866025403784439
	const taui_2 = 1.732050807568878
	for k := range l1 {
		tr2 := cc[ido-1+(3*k+1)*ido]
		tr2 = vadd(tr2, tr2)
		cr2 := vmadd(ld_ps1(taur), tr2, cc[3*k*ido])
		ch[k*ido] = vadd(cc[3*k*ido], tr2)
		ci3 := svmul(taui_2, cc[(3*k+2)*ido])
		ch[(k+l1)*ido] = vsub(cr2, ci3)
		ch[(k+2*l1)*ido] = vadd(cr2, ci3)
	}
	if ido == 1 {
		return
	}
	for k := range l1 {
		for i := 2; i < ido; i += 2 {
			ic := ido - i
			tr2 := vadd(cc[i-1+(3*k+2)*ido], cc[ic-1+(3*k+1)*ido])
			cr2 := vmadd(ld_ps1(taur), tr2, cc[i-1+3*k*ido])
			ch[i-1+k*ido] = vadd(cc[i-1+3*k*ido], tr2)
			ti2 := vsub(cc[i+(3*k+2)*ido], cc[ic+(3*k+1)*ido])
			ci2 := vmadd(ld_ps1(taur), ti2, cc[i+3*k*ido])
			ch[i+k*ido] = vadd(cc[i+3*k*ido], ti2)
			cr3 := svmul(taui, vsub(cc[i-1+(3*k+2)*ido], cc[ic-1+(3*k+1)*ido]))
			ci3 := svmul(taui, vadd(cc[i+(3*k+2)*ido], cc[ic+(3*k+1)*ido]))
			dr2 := vsub(cr2, ci3)
			dr3 := vadd(cr2, ci3)
			di2 := vadd(ci2, cr3)
			di3 := vsub(ci2, cr3)
			dr2, di2 = vcplxmul(dr2, di2, ld_ps1(wa1[i-2]), ld_ps1(wa1[i-1]))
			ch[i-1+(k+l1)*ido] = dr2
			ch[i+(k+l1)*ido] = di2
			dr3, di3 = vcplxmul(dr3, di3, ld_ps1(wa2[i-2]), ld_ps1(wa2[i-1]))
			ch[i-1+(k+2*l1)*ido] = dr3
			ch[i+(k+2*l1)*ido] = di3
		}
	}
}

//go:noinline
func radf4_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2, wa3 []float32) {
	l1ido := l1 * ido
	{
		cc_, ch_ := cc, ch
		for range l1ido {
			a0, a1 := cc[0], cc[l1ido]
			a2, a3 := cc[2*l1ido], cc[3*l1ido]
			tr1 := vadd(a1, a3)
			tr2 := vadd(a0, a2)
			ch[2*ido-1] = vsub(a0, a2)
			ch[2*ido] = vsub(a3, a1)
			ch[0] = vadd(tr1, tr2)
			ch[4*ido-1] = vsub(tr2, tr1)

			cc, ch = cc[ido:], ch[4*ido:]
		}
		cc, ch = cc_, ch_
	}
	if ido < 2 {
		return
	}
	if ido != 2 {
		for k := 0; k < l1ido; k += ido {
			pc := cc[1+k:]
			for i := 2; i < ido; i, pc = i+2, pc[2:] {
				ic := ido - i

				cr2 := pc[1*l1ido+0]
				ci2 := pc[1*l1ido+1]
				wr := ld_ps1(wa1[i-2])
				wi := ld_ps1(wa1[i-1])
				cr2, ci2 = vcplxmulconj(cr2, ci2, wr, wi)

				cr3 := pc[2*l1ido+0]
				ci3 := pc[2*l1ido+1]
				wr = ld_ps1(wa2[i-2])
				wi = ld_ps1(wa2[i-1])
				cr3, ci3 = vcplxmulconj(cr3, ci3, wr, wi)

				cr4 := pc[3*l1ido]
				ci4 := pc[3*l1ido+1]
				wr = ld_ps1(wa3[i-2])
				wi = ld_ps1(wa3[i-1])
				cr4, ci4 = vcplxmulconj(cr4, ci4, wr, wi)

				// Original comments below
				/* at this point, on sse, five of "cr2 cr3 cr4 ci2 ci3 ci4" should be loaded in registers */

				tr1 := vadd(cr2, cr4)
				tr4 := vsub(cr4, cr2)
				tr2 := vadd(pc[0], cr3)
				tr3 := vsub(pc[0], cr3)
				ch[i-1+4*k] = vadd(tr1, tr2)
				ch[ic-1+4*k+3*ido] = vsub(tr2, tr1) // at this point tr1 and tr2 can be disposed
				ti1 := vadd(ci2, ci4)
				ti4 := vsub(ci2, ci4)
				ch[i-1+4*k+2*ido] = vadd(ti4, tr3)
				ch[ic-1+4*k+1*ido] = vsub(tr3, ti4) // dispose tr3, ti4
				ti2 := vadd(pc[1], ci3)
				ti3 := vsub(pc[1], ci3)
				ch[i+4*k] = vadd(ti1, ti2)
				ch[ic+4*k+3*ido] = vsub(ti1, ti2)
				ch[i+4*k+2*ido] = vadd(tr4, ti3)
				ch[ic+4*k+1*ido] = vsub(tr4, ti3)
			}
		}
		if ido%2 == 1 {
			return
		}
	}
	const nrsqrt2 = -0.7071067811865475
	for k := 0; k < l1ido; k += ido {
		a, b := cc[ido-1+k+3*l1ido], cc[ido-1+k+l1ido]
		c, d := cc[ido-1+k+2*l1ido], cc[ido-1+k]
		ti1 := svmul(nrsqrt2, vadd(a, b))
		tr1 := svmul(nrsqrt2, vsub(b, a))
		ch[ido-1+4*k] = vadd(tr1, c)
		ch[ido-1+4*k+2*ido] = vsub(c, tr1)
		ch[4*k+1*ido] = vsub(ti1, d)
		ch[4*k+3*ido] = vadd(ti1, d)
	}
}

//go:noinline
func radb4_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2, wa3 []float32) {
	l1ido := l1 * ido
	{
		cc_, ch_ := cc, ch
		for range l1ido {
			a, b := cc[0], cc[4*l1ido-1]
			c, d := cc[2*l1ido], cc[2*l1ido-1]
			tr3 := svmul(2, d)
			tr2 := vadd(a, b)
			tr1 := vsub(a, b)
			tr4 := svmul(2, c)
			ch[0*l1ido] = vadd(tr2, tr3)
			ch[2*l1ido] = vsub(tr2, tr3)
			ch[1*l1ido] = vsub(tr1, tr4)
			ch[3*l1ido] = vadd(tr1, tr4)

			cc, ch = cc[ido:], ch[4*ido:]
		}
		cc, ch = cc_, ch_
	}
	if ido != 2 {
		for k := 0; k < l1ido; k += ido {
			pc := cc[1+4*k:]
			ph := cc[1+k:]
			for i := 2; i < ido; i += 2 {
				tr1 := vsub(pc[i], pc[4*ido-i])
				tr2 := vadd(pc[i], pc[4*ido-i])
				ti4 := vsub(pc[2*ido+i], pc[2*ido-i])
				tr3 := vadd(pc[2*ido+i], pc[2*ido-i])
				ph[0] = vadd(tr2, tr3)
				cr3 := vsub(tr2, tr3)

				ti3 := vsub(pc[2*ido+i+1], pc[2*ido-i+1])
				tr4 := vadd(pc[2*ido+i+1], pc[2*ido-i+1])
				cr2 := vsub(tr1, tr4)
				cr4 := vadd(tr1, tr4)

				ti1 := vadd(pc[i+1], pc[4*ido-i+1])
				ti2 := vsub(pc[i+1], pc[4*ido-i+1])

				oph := ph
				ph[1] = vadd(ti2, ti3)
				ph = ph[l1ido:]
				ci3 := vsub(ti2, ti3)
				ci2 := vadd(ti1, ti4)
				ci4 := vsub(ti1, ti4)
				cr2, ci2 = vcplxmul(cr2, ci2, ld_ps1(wa1[i-2]), ld_ps1(wa1[i-1]))
				ph[0] = cr2
				ph[1] = ci2
				ph = ph[l1ido:]
				cr3, ci3 = vcplxmul(cr3, ci3, ld_ps1(wa2[i-2]), ld_ps1(wa2[i-1]))
				ph[0] = cr3
				ph[1] = ci3
				ph = ph[l1ido:]
				cr4, ci4 = vcplxmul(cr4, ci4, ld_ps1(wa3[i-2]), ld_ps1(wa3[i-1]))
				ph[0] = cr4
				ph[1] = ci4
				ph = oph[2:]
			}
		}
		if ido%2 == 1 {
			return
		}
	}
	const nsqrt2 = -1.414213562373095
	for k := 0; k < l1ido; k += ido {
		i0 := 4*k + ido
		c, d := cc[i0+2*ido-1], cc[i0-1]
		a, b := cc[i0+2*ido+0], cc[i0+0]
		tr1 := vsub(c, d)
		tr2 := vadd(c, d)
		ti1 := vadd(b, a)
		ti2 := vsub(b, a)
		ch[ido-1+k+0*l1ido] = vadd(tr2, tr2)
		ch[ido-1+k+1*l1ido] = svmul(-nsqrt2, vsub(ti1, tr1))
		ch[ido-1+k+2*l1ido] = vadd(ti2, ti2)
		ch[ido-1+k+3*l1ido] = svmul(-nsqrt2, vadd(ti1, tr1))
	}
}

//go:noinline
func radf5_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2, wa3, wa4 []float32) {
	const tr11 = .309016994374947
	const ti11 = .951056516295154
	const tr12 = -.809016994374947
	const ti12 = .587785252292473

	ch_offset := 1 + ido*6
	cc_offset := 1 + ido*(1+l1)

	cc_ref := func(a_1, a_2, a_3 int) simd.Float32s { return cc[((a_3)*l1+(a_2))*ido+a_1-cc_offset] }
	ch_ref := func(a_1, a_2, a_3 int) *simd.Float32s { return &ch[((a_3)*5+(a_2))*ido+a_1-ch_offset] }

	for k := 1; k <= l1; k++ {
		cr2 := vadd(cc_ref(1, k, 5), cc_ref(1, k, 2))
		ci5 := vsub(cc_ref(1, k, 5), cc_ref(1, k, 2))
		cr3 := vadd(cc_ref(1, k, 4), cc_ref(1, k, 3))
		ci4 := vsub(cc_ref(1, k, 4), cc_ref(1, k, 3))
		*ch_ref(1, 1, k) = vadd(cc_ref(1, k, 1), vadd(cr2, cr3))
		*ch_ref(ido, 2, k) = vadd(cc_ref(1, k, 1), vadd(svmul(tr11, cr2), svmul(tr12, cr3)))
		*ch_ref(1, 3, k) = vadd(svmul(ti11, ci5), svmul(ti12, ci4))
		*ch_ref(ido, 4, k) = vadd(cc_ref(1, k, 1), vadd(svmul(tr12, cr2), svmul(tr11, cr3)))
		*ch_ref(1, 5, k) = vsub(svmul(ti12, ci5), svmul(ti11, ci4))
		//printf("pffft: radf5, k=%d ch_ref=%f, ci4=%f\n", k, ch_ref(1, 5, k), ci4);
	}
	if ido == 1 {
		return
	}
	idp2 := ido + 2
	for k := 1; k <= l1; k++ {
		for i := 3; i <= ido; i += 2 {
			ic := idp2 - i
			dr2 := ld_ps1(wa1[i-3])
			di2 := ld_ps1(wa1[i-2])
			dr3 := ld_ps1(wa2[i-3])
			di3 := ld_ps1(wa2[i-2])
			dr4 := ld_ps1(wa3[i-3])
			di4 := ld_ps1(wa3[i-2])
			dr5 := ld_ps1(wa4[i-3])
			di5 := ld_ps1(wa4[i-2])
			dr2, di2 = vcplxmulconj(dr2, di2, cc_ref(i-1, k, 2), cc_ref(i, k, 2))
			dr3, di3 = vcplxmulconj(dr3, di3, cc_ref(i-1, k, 3), cc_ref(i, k, 3))
			dr4, di4 = vcplxmulconj(dr4, di4, cc_ref(i-1, k, 4), cc_ref(i, k, 4))
			dr5, di5 = vcplxmulconj(dr5, di5, cc_ref(i-1, k, 5), cc_ref(i, k, 5))
			cr2 := vadd(dr2, dr5)
			ci5 := vsub(dr5, dr2)
			cr5 := vsub(di2, di5)
			ci2 := vadd(di2, di5)
			cr3 := vadd(dr3, dr4)
			ci4 := vsub(dr4, dr3)
			cr4 := vsub(di3, di4)
			ci3 := vadd(di3, di4)
			*ch_ref(i-1, 1, k) = vadd(cc_ref(i-1, k, 1), vadd(cr2, cr3))
			*ch_ref(i, 1, k) = vsub(cc_ref(i, k, 1), vadd(ci2, ci3))
			tr2 := vadd(cc_ref(i-1, k, 1), vadd(svmul(tr11, cr2), svmul(tr12, cr3)))
			ti2 := vsub(cc_ref(i, k, 1), vadd(svmul(tr11, ci2), svmul(tr12, ci3)))
			tr3 := vadd(cc_ref(i-1, k, 1), vadd(svmul(tr12, cr2), svmul(tr11, cr3)))
			ti3 := vsub(cc_ref(i, k, 1), vadd(svmul(tr12, ci2), svmul(tr11, ci3)))
			tr5 := vadd(svmul(ti11, cr5), svmul(ti12, cr4))
			ti5 := vadd(svmul(ti11, ci5), svmul(ti12, ci4))
			tr4 := vsub(svmul(ti12, cr5), svmul(ti11, cr4))
			ti4 := vsub(svmul(ti12, ci5), svmul(ti11, ci4))
			*ch_ref(i-1, 3, k) = vsub(tr2, tr5)
			*ch_ref(ic-1, 2, k) = vadd(tr2, tr5)
			*ch_ref(i, 3, k) = vadd(ti2, ti5)
			*ch_ref(ic, 2, k) = vsub(ti5, ti2)
			*ch_ref(i-1, 5, k) = vsub(tr3, tr4)
			*ch_ref(ic-1, 4, k) = vadd(tr3, tr4)
			*ch_ref(i, 5, k) = vadd(ti3, ti4)
			*ch_ref(ic, 4, k) = vsub(ti4, ti3)
		}
	}
}

//go:noinline
func radb5_ps(ido, l1 int, cc, ch []simd.Float32s, wa1, wa2, wa3, wa4 []float32) {
	const tr11 = .309016994374947
	const ti11 = .951056516295154
	const tr12 = -.809016994374947
	const ti12 = .587785252292473

	ch_offset := 1 + ido*(1+l1)
	cc_offset := 1 + ido*6

	cc_ref := func(a_1, a_2, a_3 int) simd.Float32s { return cc[((a_3)*5+(a_2))*ido+a_1-cc_offset] }
	ch_ref := func(a_1, a_2, a_3 int) *simd.Float32s { return &ch[((a_3)*l1+(a_2))*ido+a_1-ch_offset] }

	for k := 1; k <= l1; k++ {
		ti5 := vadd(cc_ref(1, 3, k), cc_ref(1, 3, k))
		ti4 := vadd(cc_ref(1, 5, k), cc_ref(1, 5, k))
		tr2 := vadd(cc_ref(ido, 2, k), cc_ref(ido, 2, k))
		tr3 := vadd(cc_ref(ido, 4, k), cc_ref(ido, 4, k))
		*ch_ref(1, k, 1) = vadd(cc_ref(1, 1, k), vadd(tr2, tr3))
		cr2 := vadd(cc_ref(1, 1, k), vadd(svmul(tr11, tr2), svmul(tr12, tr3)))
		cr3 := vadd(cc_ref(1, 1, k), vadd(svmul(tr12, tr2), svmul(tr11, tr3)))
		ci5 := vadd(svmul(ti11, ti5), svmul(ti12, ti4))
		ci4 := vsub(svmul(ti12, ti5), svmul(ti11, ti4))
		*ch_ref(1, k, 2) = vsub(cr2, ci5)
		*ch_ref(1, k, 3) = vsub(cr3, ci4)
		*ch_ref(1, k, 4) = vadd(cr3, ci4)
		*ch_ref(1, k, 5) = vadd(cr2, ci5)
	}
	if ido == 1 {
		return
	}
	idp2 := ido + 2
	for k := 1; k <= l1; k++ {
		for i := 3; i <= ido; i += 2 {
			ic := idp2 - i
			ti5 := vadd(cc_ref(i, 3, k), cc_ref(ic, 2, k))
			ti2 := vsub(cc_ref(i, 3, k), cc_ref(ic, 2, k))
			ti4 := vadd(cc_ref(i, 5, k), cc_ref(ic, 4, k))
			ti3 := vsub(cc_ref(i, 5, k), cc_ref(ic, 4, k))
			tr5 := vsub(cc_ref(i-1, 3, k), cc_ref(ic-1, 2, k))
			tr2 := vadd(cc_ref(i-1, 3, k), cc_ref(ic-1, 2, k))
			tr4 := vsub(cc_ref(i-1, 5, k), cc_ref(ic-1, 4, k))
			tr3 := vadd(cc_ref(i-1, 5, k), cc_ref(ic-1, 4, k))
			*ch_ref(i-1, k, 1) = vadd(cc_ref(i-1, 1, k), vadd(tr2, tr3))
			*ch_ref(i, k, 1) = vadd(cc_ref(i, 1, k), vadd(ti2, ti3))
			cr2 := vadd(cc_ref(i-1, 1, k), vadd(svmul(tr11, tr2), svmul(tr12, tr3)))
			ci2 := vadd(cc_ref(i, 1, k), vadd(svmul(tr11, ti2), svmul(tr12, ti3)))
			cr3 := vadd(cc_ref(i-1, 1, k), vadd(svmul(tr12, tr2), svmul(tr11, tr3)))
			ci3 := vadd(cc_ref(i, 1, k), vadd(svmul(tr12, ti2), svmul(tr11, ti3)))
			cr5 := vadd(svmul(ti11, tr5), svmul(ti12, tr4))
			ci5 := vadd(svmul(ti11, ti5), svmul(ti12, ti4))
			cr4 := vsub(svmul(ti12, tr5), svmul(ti11, tr4))
			ci4 := vsub(svmul(ti12, ti5), svmul(ti11, ti4))
			dr3 := vsub(cr3, ci4)
			dr4 := vadd(cr3, ci4)
			di3 := vadd(ci3, cr4)
			di4 := vsub(ci3, cr4)
			dr5 := vadd(cr2, ci5)
			dr2 := vsub(cr2, ci5)
			di5 := vsub(ci2, cr5)
			di2 := vadd(ci2, cr5)
			dr2, di2 = vcplxmul(dr2, di2, ld_ps1(wa1[i-3]), ld_ps1(wa1[i-2]))
			dr3, di3 = vcplxmul(dr3, di3, ld_ps1(wa2[i-3]), ld_ps1(wa2[i-2]))
			dr4, di4 = vcplxmul(dr4, di4, ld_ps1(wa3[i-3]), ld_ps1(wa3[i-2]))
			dr5, di5 = vcplxmul(dr5, di5, ld_ps1(wa4[i-3]), ld_ps1(wa4[i-2]))

			*ch_ref(i-1, k, 2) = dr2
			*ch_ref(i, k, 2) = di2
			*ch_ref(i-1, k, 3) = dr3
			*ch_ref(i, k, 3) = di3
			*ch_ref(i-1, k, 4) = dr4
			*ch_ref(i, k, 4) = di4
			*ch_ref(i-1, k, 5) = dr5
			*ch_ref(i, k, 5) = di5
		}
	}
}
