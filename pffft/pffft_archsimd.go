//go:build goexperiment.simd

package pffft

import "simd"

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
				vcplxmulconj(tr2, ti2, ld_ps1(wa1[i-2]), ld_ps1(wa1[i-1]))
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
