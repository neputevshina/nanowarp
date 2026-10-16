//go:build goexperiment.simd

package gopffft

import (
	"fmt"
	"math"
	"simd"
	"unsafe"
)

func init() {
	if simd.VectorBitSize() != 128 {
		panic(`gopffft currently supports only 128-bit vectors, run the program with GODEBUG=simd=+128`)
	}
}

type reg = simd.Float32s

func getAsAligned(in []float32) []reg {
	nvec := simd.VectorBitSize() / 32
	if len(in)%nvec != 0 {
		panic(`gopffft: input slice length must be divided by platform's vector size`)
	}
	c := uintptr(simd.VectorBitSize() / 8)
	if uintptr(unsafe.Pointer(unsafe.SliceData(in)))&(c-1) != 0 {
		panic(`gopffft: unaligned input`)
	}
	return unsafe.Slice((*reg)(unsafe.Pointer(unsafe.SliceData(in))), len(in)/nvec)
}

func degrade(in []reg) []float32 {
	nvec := simd.VectorBitSize() / 32
	return unsafe.Slice((*float32)(unsafe.Pointer(unsafe.SliceData(in))), len(in)*nvec)
}

func align(in []float32) []float32 {
	nvec := simd.VectorBitSize() / 32
	rem := uintptr(unsafe.Pointer(unsafe.SliceData(in))) % uintptr(nvec)
	in = in[nvec-int(rem):]
	_ = getAsAligned(in)
	return in
}

// v4sf tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br);
// ar=VSUB(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VADD(ai,tmp); }
func vcplxmul(ar, ai, br, bi reg) (reg, reg) {
	ar = ar.Mul(br).Sub(ai.Mul(bi))
	ai = ai.Mul(br).Add(ar.Mul(bi))
	return ar, ai
}

// v4sf tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br);
// ar=VADD(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VSUB(ai,tmp);
func vcplxmulconj(ar, ai, br, bi reg) (reg, reg) {
	ar = ar.Mul(br).Add(ai.Mul(bi))
	ai = ai.Mul(br).Sub(ar.Mul(bi))
	return ar, ai
}

func svmul(s float32, ps reg) reg {
	return simd.BroadcastFloat32s(s).Mul(ps)
}

func vadd(a, b reg) reg {
	return a.Add(b)
}

func vsub(a, b reg) reg {
	return a.Sub(b)
}

func vmul(a, b reg) reg {
	return a.Mul(b)
}

func vmadd(a, b, c reg) reg {
	return a.MulAdd(b, c)
}

func ld_ps1(s float32) reg {
	return simd.BroadcastFloat32s(s)
}

// TODO Shuffle is not yet implemented in simd and in archsimd it's only for words and
// only for AVX.
func interleave2(in1, in2 reg) (out1, out2 reg) {
	var x, y, out [32]float32
	in1.Store(x[:])
	in2.Store(y[:])
	out[0] = x[0]
	out[1] = y[0]
	out[2] = x[1]
	out[3] = y[1]
	out[4] = x[2]
	out[5] = y[2]
	out[6] = x[3]
	out[7] = y[3]
	out1 = simd.LoadFloat32s(out[:4])
	out2 = simd.LoadFloat32s(out[4:8])
	return
}

// TODO
func uninterleave2(in1, in2 reg) (out1, out2 reg) {
	var x, y, out [32]float32
	in1.Store(out[:])
	in2.Store(out[4:])
	x[0] = out[0]
	y[0] = out[1]
	x[1] = out[2]
	y[1] = out[3]
	x[2] = out[4]
	y[2] = out[5]
	x[3] = out[6]
	y[3] = out[7]
	out1 = simd.LoadFloat32s(x[:])
	out2 = simd.LoadFloat32s(y[:])
	return
}

// TODO
func vtranspose4(a, b, c, d reg) (a1, b1, c1, d1 reg) {
	var x, y [32]float32
	a.Store(x[:])
	b.Store(x[4:])
	c.Store(x[8:])
	d.Store(x[12:])
	for i := range 4 {
		y[i*4+0] = x[i+0*4]
		y[i*4+1] = x[i+1*4]
		y[i*4+2] = x[i+2*4]
		y[i*4+3] = x[i+3*4]
	}
	a1.Store(y[:])
	b1.Store(y[4:])
	c1.Store(y[8:])
	d1.Store(y[12:])
	return
}

// TODO
func vswaphl(a, b reg) (out reg) {
	var x, y [32]float32
	a.Store(x[:])
	b.Store(x[4:])
	y[0] = x[0]
	y[1] = x[1]
	y[2] = x[4]
	y[3] = x[5]
	out.Store(y[:])
	return
}

//go:noinline
func passf2_ps(ido, l1 int, cc, ch []reg, wa1 []float32, fsign float32) {
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
func passf3_ps(ido, l1 int, cc, ch []reg, wa1, wa2 []float32, fsign float32) {
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
func passf4_ps(ido, l1 int, cc, ch []reg, wa1, wa2, wa3 []float32, fsign float32) {
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
func passf5_ps(ido, l1 int, cc, ch []reg, wa1, wa2, wa3, wa4 []float32, fsign float32) {
	const tr11 = .309016994374947
	ti11 := .951056516295154 * fsign
	const tr12 = -.809016994374947
	ti12 := .587785252292473 * fsign

	if ido <= 2 {
		panic(`fatal: passf5_ps: ido <= 2`)
	}
	cc_ref := func(a_1, a_2 int) reg { return cc[(a_2-1)*ido+a_1+1] }
	ch_ref := func(a_1, a_3 int) *reg { return &ch[(a_3-1)*l1*ido+a_1+1] }

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
func radf2_ps(ido, l1 int, cc, ch []reg, wa1 []float32) {
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
func radb2_ps(ido, l1 int, cc, ch []reg, wa1 []float32) {
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
func radf3_ps(ido, l1 int, cc, ch []reg, wa1, wa2 []float32) {
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
func radb3_ps(ido, l1 int, cc, ch []reg, wa1, wa2 []float32) {
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
func radf4_ps(ido, l1 int, cc, ch []reg, wa1, wa2, wa3 []float32) {
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
func radb4_ps(ido, l1 int, cc, ch []reg, wa1, wa2, wa3 []float32) {
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
func radf5_ps(ido, l1 int, cc, ch []reg, wa1, wa2, wa3, wa4 []float32) {
	const tr11 = .309016994374947
	const ti11 = .951056516295154
	const tr12 = -.809016994374947
	const ti12 = .587785252292473

	ch_offset := 1 + ido*6
	cc_offset := 1 + ido*(1+l1)

	cc_ref := func(a_1, a_2, a_3 int) reg { return cc[((a_3)*l1+(a_2))*ido+a_1-cc_offset] }
	ch_ref := func(a_1, a_2, a_3 int) *reg { return &ch[((a_3)*5+(a_2))*ido+a_1-ch_offset] }

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
func radb5_ps(ido, l1 int, cc, ch []reg, wa1, wa2, wa3, wa4 []float32) {
	const tr11 = .309016994374947
	const ti11 = .951056516295154
	const tr12 = -.809016994374947
	const ti12 = .587785252292473

	ch_offset := 1 + ido*(1+l1)
	cc_offset := 1 + ido*6

	cc_ref := func(a_1, a_2, a_3 int) reg { return cc[((a_3)*5+(a_2))*ido+a_1-cc_offset] }
	ch_ref := func(a_1, a_2, a_3 int) *reg { return &ch[((a_3)*l1+(a_2))*ido+a_1-ch_offset] }

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

//go:noinline
func rfftf1_ps(n int, input_readonly, work1, work2 []simd.Float32s, wa []float32, ifac []int) []simd.Float32s {
	in := input_readonly
	out := work2
	if &in[0] == &work2[0] {
		out = work1
	}
	nf := ifac[1]
	l2 := n
	iw := n - 1
	if &in[0] == &out[0] || &work1[0] == &work2[0] {
		panic(`fatal: in-place work or out`)
	}
	for k1 := 1; k1 <= nf; k1++ {
		kh := nf - k1
		ip := ifac[kh+2]
		l1 := l2 / ip
		ido := n / l2
		iw -= (ip - 1) * ido
		switch ip {
		case 5:
			ix2 := iw + ido
			ix3 := ix2 + ido
			ix4 := ix3 + ido
			radf5_ps(ido, l1, in, out, wa[iw:], wa[ix2:], wa[ix3:], wa[ix4:])
		case 4:
			ix2 := iw + ido
			ix3 := ix2 + ido
			radf4_ps(ido, l1, in, out, wa[iw:], wa[ix2:], wa[ix3:])
		case 3:
			ix2 := iw + ido
			radf3_ps(ido, l1, in, out, wa[iw:], wa[ix2:])
		case 2:
			radf2_ps(ido, l1, in, out, wa[iw:])
		default:
			panic(`fatal: wrong radix (ip)`)
		}
		l2 = l1
		if &out[0] == &work2[0] {
			out, in = work1, work2
		} else {
			out, in = work2, work1
		}
	}
	if &in[0] == &input_readonly[0] {
		panic(`fatal: dirty input`)
	}
	return in
}

//go:noinline
func rfftb1_ps(n int, input_readonly, work1, work2 []simd.Float32s, wa []float32, ifac []int) []simd.Float32s {
	in := input_readonly
	out := work2
	if &in[0] == &work2[0] {
		out = work1
	}
	nf := ifac[1]
	l1 := 1
	iw := 0
	if &in[0] == &out[0] {
		panic(`fatal: in-place out`)
	}
	for k1 := 1; k1 <= nf; k1++ {
		ip := ifac[k1+1]
		l2 := ip * l1
		ido := n / l2
		iw -= (ip - 1) * ido
		switch ip {
		case 5:
			ix2 := iw + ido
			ix3 := ix2 + ido
			ix4 := ix3 + ido
			radb5_ps(ido, l1, in, out, wa[iw:], wa[ix2:], wa[ix3:], wa[ix4:])
		case 4:
			ix2 := iw + ido
			ix3 := ix2 + ido
			radb4_ps(ido, l1, in, out, wa[iw:], wa[ix2:], wa[ix3:])
		case 3:
			ix2 := iw + ido
			radb3_ps(ido, l1, in, out, wa[iw:], wa[ix2:])
		case 2:
			radb2_ps(ido, l1, in, out, wa[iw:])
		default:
			panic(`fatal: wrong radix (ip)`)
		}
		l1 = l2
		iw += (ip - 1) * ido

		if &out[0] == &work2[0] {
			out, in = work1, work2
		} else {
			out, in = work2, work1
		}
	}
	if &in[0] == &input_readonly[0] {
		panic(`fatal: dirty input`)
	}
	return in
}

/* max number of integer factors for the decomposition, +2 */
const IfacMaxSize = 25

func decompose(n int, ifac, ntryh []int) int {
	nl := n
	nf := 0
	for j := 0; ntryh[j] != 0; j++ {
		ntry := ntryh[j]
		for nl != 1 {
			nq := nl / ntry
			nr := nl - ntry*nq
			if nr == 0 {
				if !(2+nf < IfacMaxSize) {
					panic(`fatal: too composite n`)
				}
				ifac[2+nf] = ntry
				nf++
				nl = nq
				if ntry == 2 && nf != 1 {
					for i := 2; i <= nf; i++ {
						ib := nf - i + 2
						ifac[ib+1] = ifac[ib]
					}
					ifac[2] = 2
				}
			} else {
				break
			}
		}
	}
	ifac[0] = n
	ifac[1] = nf
	return nf
}

func rffti1_ps(n int, wa []float32, ifac []int) {
	ntryh := []int{4, 2, 3, 5, 0}

	nf := decompose(n, ifac, ntryh)
	argh := 2 * math.Pi / float32(n)
	is := 0
	l1 := 1
	for k1 := 1; k1 <= nf-1; k1++ {
		ip := ifac[k1+1]
		ld := 0
		l2 := l1 * ip
		ido := n / l2
		for range ip - 1 {
			i := is
			fi := 0
			ld += l1
			argld := float32(ld) * argh
			for ii := 3; i <= ido; ii += 2 {
				i += 2
				fi += 1
				a := float64(float32(fi) * argld)
				wa[i-2] = float32(math.Cos(a))
				wa[i-1] = float32(math.Sin(a))
			}
			is += ido
		}
		l1 = l2
	}
}

func cffti1_ps(n int, wa []float32, ifac []int) {
	ntryh := []int{5, 3, 4, 2, 0}

	nf := decompose(n, ifac, ntryh)
	argh := 2 * math.Pi / float32(n)
	i := 1
	l1 := 1
	for k1 := 1; k1 <= nf; k1++ {
		ip := ifac[k1+1]
		ld := 0
		l2 := l1 * ip
		ido := n / l2
		idot := ido + ido + 2
		for range ip - 1 {
			i1 := i
			fi := 0
			wa[i-1] = 1
			wa[i] = 0
			ld += l1
			argld := float32(ld) * argh
			for ii := 4; i <= idot; ii += 2 {
				i += 2
				fi += 1
				a := float64(float32(fi) * argld)
				wa[i-1] = float32(math.Cos(a))
				wa[i] = float32(math.Sin(a))
			}
			if ip > 5 {
				wa[i1-i] = wa[i-1]
				wa[i1] = wa[i]
			}
		}
		l1 = l2
	}
}

func cfftf1_ps(n int, input_readonly, work1, work2 []simd.Float32s, wa []float32, ifac []int, isign int) []simd.Float32s {
	in := input_readonly
	out := work2
	if &in[0] == &work2[0] {
		out = work1
	}
	nf := ifac[1]
	l1 := 1
	iw := 0
	if &in[0] == &out[0] || &work1[0] == &work2[0] {
		panic(`fatal: in-place work or out`)
	}
	for k1 := 2; k1 <= nf+1; k1++ {
		ip := ifac[k1]
		l2 := ip * l1
		ido := n / l2
		idot := ido + ido
		switch ip {
		case 5:
			ix2 := iw + idot
			ix3 := ix2 + idot
			ix4 := ix3 + idot
			passf5_ps(ido, l1, in, out, wa[iw:], wa[ix2:], wa[ix3:], wa[ix4:], float32(isign))
		case 4:
			ix2 := iw + idot
			ix3 := ix2 + idot
			passf4_ps(ido, l1, in, out, wa[iw:], wa[ix2:], wa[ix3:], float32(isign))
		case 3:
			ix2 := iw + idot
			passf3_ps(ido, l1, in, out, wa[iw:], wa[ix2:], float32(isign))
		case 2:
			passf2_ps(ido, l1, in, out, wa[iw:], float32(isign))
		default:
			panic(`fatal: wrong radix (ip)`)
		}
		l1 = l2
		iw += (ip - 1) * idot

		if &out[0] == &work2[0] {
			out, in = work1, work2
		} else {
			out, in = work2, work1
		}
	}
	if &in[0] == &input_readonly[0] {
		panic(`fatal: dirty input`)
	}
	return in
}

func pffftNewSetup(n int, transform transform) (s *pffftSetup) {
	SIMD_SZ := simd.VectorBitSize() / 32
	// validate N for negative values or potential int overflow
	if n < 0 {
		return nil
	}
	if n > 1<<26 {
		panic(`pffftNewSetup: nfft too large`)
	}
	s = &pffftSetup{}
	/* unfortunately, the fft size must be a multiple of 16 for complex FFTs
	   and 32 for real FFTs -- a lot of stuff would need to be rewritten to
	   handle other cases (or maybe just switch to a scalar fft, I don't know..) */
	if transform == transformReal && n%(2*SIMD_SZ*SIMD_SZ) != 0 {
		panic(fmt.Sprint(`pffftNewSetup: nfft must be a multiple of`, 2*SIMD_SZ*SIMD_SZ))
	}
	if transform == transformComplex && n%(SIMD_SZ*SIMD_SZ) != 0 {
		panic(fmt.Sprint(`pffftNewSetup: nfft must be a multiple of`, SIMD_SZ*SIMD_SZ))
	}
	s.n = n
	s.transform = transform
	/* nb of complex simd vectors */
	if transform == transformReal {
		s.ncvec = n / 2 / SIMD_SZ
	} else {
		s.ncvec = n / SIMD_SZ
	}
	s.data = make([]reg, 2*s.ncvec)
	s.e = s.data
	// FIXME A hack, the real way is to iterate over s.e taking bites of size SIMD_SZ.
	// Or simply to never use vectors directly, since this is technically a wrong usage.
	e := unsafe.Slice((*float32)(unsafe.Pointer(unsafe.SliceData(s.e))), len(s.data)*SIMD_SZ)
	s.twiddle = make([]float32, 2*n)

	for k := range s.ncvec {
		i, j := k/SIMD_SZ, k%SIMD_SZ
		for m := range SIMD_SZ - 1 {
			A := -2 * math.Pi * float64((m+1)*k) / float64(n)
			// TODO Will this even work when SIMD_SZ is not 4 or 1?
			e[(2*(i*3+m)+0)*SIMD_SZ+j] = float32(math.Cos(A))
			e[(2*(i*3+m)+1)*SIMD_SZ+j] = float32(math.Sin(A))
		}
	}

	if transform == transformReal {
		rffti1_ps(n/SIMD_SZ, s.twiddle, s.ifac[:])
	} else {
		cffti1_ps(n/SIMD_SZ, s.twiddle, s.ifac[:])
	}

	/* check that N is decomposable with allowed prime factors */
	m := 1
	for k := range s.ifac[1] {
		m *= s.ifac[2*k]
	}
	if m != n/SIMD_SZ {
		s = nil
	}

	return
}

func reversed_copy(n int, in []reg, stride int, out []reg) {
	forw := func() (a0, a1 reg) {
		a0, a1 = interleave2(in[0], in[1])
		in = in[stride:]
		return
	}
	bakw := func(a0, a1 reg) {
		out[len(out)-1] = vswaphl(a0, a1)
		out = out[:len(out)-1]
	}

	g0, g1 := forw()
	bakw(g0, g1)
	for range n - 1 {
		h0, h1 := forw()
		bakw(g1, h0)
		bakw(h0, h1)
		g1 = h1
	}
	bakw(g1, g0)
}

func unreversed_copy(n int, in []reg, out []reg, stride int) {
	g0, g1 := in[0], in[0]
	in = in[1:]

	adv := func(h0, h1 reg) reg {
		g1 = vswaphl(g1, h0)
		h0 = vswaphl(h0, h1)
		out[0], out[1] = uninterleave2(h0, g1)
		return h1
	}

	for range n - 1 {
		g1 = adv(in[0], in[1])
		in = in[2:]
		out = out[stride:]
	}
	_ = adv(in[0], g0)
}

func pffft_zreorder(s *pffftSetup, in, out []float32, forward bool) {
	SIMD_SZ := simd.VectorBitSize() / 32
	n, ncvec := s.n, s.ncvec
	vin, vout := getAsAligned(in), getAsAligned(out)
	if &in[0] == &out[0] {
		panic(`gopffft: pffft_zreorder: in and out can't alias`)
	}
	dk := n / 32

	switch {
	case s.transform == transformReal && forward:
		for k := range dk {
			vout[2*(0*dk+k)+0], vout[2*(0*dk+k)+1] = interleave2(vin[k*8+0], vin[k*8+1])
			vout[2*(2*dk+k)+0], vout[2*(2*dk+k)+1] = interleave2(vin[k*8+4], vin[k*8+5])
		}
		// NOTE(neputevshina): No cast, reworked to point directly to vin and vout.
		// n is always a multiple of SIMD_SZ.
		reversed_copy(dk, vin[2:], 8, vout[n/2/SIMD_SZ:])
		reversed_copy(dk, vin[6:], 8, vout[n/SIMD_SZ:])
	case s.transform == transformReal && !forward:
		for k := range dk {
			vout[k*8+0], vout[k*8+1] = uninterleave2(vin[2*(0*dk+k)+0], vin[2*(0*dk+k)+1])
			vout[k*8+4], vout[k*8+5] = uninterleave2(vin[2*(2*dk+k)+0], vin[2*(2*dk+k)+1])
		}
		// NOTE(neputevshina): No cast, reworked to point directly to vin and vout.
		// n is always a multiple of SIMD_SZ.
		unreversed_copy(dk, vin[n/4/SIMD_SZ:], vout[n/SIMD_SZ-6:], -8)
		unreversed_copy(dk, vin[3*n/4/SIMD_SZ:], vout[n/SIMD_SZ-2:], -8)
	case s.transform == transformComplex && forward:
		for k := range ncvec {
			kk := (k / 4) + (k%4)*(ncvec/4)
			vout[kk*2], vout[kk*2+1] = interleave2(vin[k*2], vin[k*2+1])
		}
	case s.transform == transformComplex && !forward:
		for k := range ncvec {
			kk := (k / 4) + (k%4)*(ncvec/4)
			vout[kk*2], vout[kk*2+1] = uninterleave2(vin[k*2], vin[k*2+1])
		}
	default:
		panic(`fatal: unknown s.transform`)
	}
}

func pffft_cplx_finalize(ncvec int, in, out, e []reg) {
	SIMD_SZ := simd.VectorBitSize() / 32
	dk := ncvec / SIMD_SZ // number of 4x4 matrix blocks
	if &in[0] == &out[0] {
		panic(`gopffft: pffft_cplx_finalize: in and out can't alias`)
	}
	for k := range dk {
		r0, i0 := in[8*k+0], in[8*k+1]
		r1, i1 := in[8*k+2], in[8*k+3]
		r2, i2 := in[8*k+4], in[8*k+5]
		r3, i3 := in[8*k+6], in[8*k+7]
		r0, r1, r2, r3 = vtranspose4(r0, r1, r2, r3)
		i0, i1, i2, i3 = vtranspose4(i0, i1, i2, i3)
		r1, i1 = vcplxmul(r1, i1, e[k*6+0], e[k*6+1])
		r2, i2 = vcplxmul(r2, i2, e[k*6+2], e[k*6+3])
		r3, i3 = vcplxmul(r3, i3, e[k*6+4], e[k*6+5])

		sr0, dr0 := vadd(r0, r2), vsub(r0, r2)
		sr1, dr1 := vadd(r1, r3), vsub(r1, r3)
		si0, di0 := vadd(i0, i2), vsub(i0, i2)
		si1, di1 := vadd(i1, i3), vsub(i1, i3)

		/*
		   transformation for each column is:

		   [1   1   1   1   0   0   0   0]   [r0]
		   [1   0  -1   0   0  -1   0   1]   [r1]
		   [1  -1   1  -1   0   0   0   0]   [r2]
		   [1   0  -1   0   0   1   0  -1]   [r3]
		   [0   0   0   0   1   1   1   1] * [i0]
		   [0   1   0  -1   1   0  -1   0]   [i1]
		   [0   0   0   0   1  -1   1  -1]   [i2]
		   [0  -1   0   1   1   0  -1   0]   [i3]
		*/

		r0, i0 = vadd(sr0, sr1), vadd(si0, si1)
		r1, i1 = vadd(dr0, di1), vsub(di0, dr1)
		r2, i2 = vsub(sr0, sr1), vsub(si0, si1)
		r3, i3 = vsub(dr0, di1), vadd(di0, dr1)

		out[0], out[1], out[2], out[3] = r0, i0, r1, i1
		out[4], out[5], out[6], out[7] = r2, i2, r3, i3
		out = out[8:]
	}
}

func pffft_cplx_preprocess(ncvec int, in, out, e []reg) {
	SIMD_SZ := simd.VectorBitSize() / 32
	dk := ncvec / SIMD_SZ // number of 4x4 matrix blocks
	if &in[0] == &out[0] {
		panic(`gopffft: pffft_cplx_preprocess: in and out can't alias`)
	}
	for k := range dk {
		r0, i0 := in[8*k+0], in[8*k+1]
		r1, i1 := in[8*k+2], in[8*k+3]
		r2, i2 := in[8*k+4], in[8*k+5]
		r3, i3 := in[8*k+6], in[8*k+7]

		sr0, dr0 := vadd(r0, r2), vsub(r0, r2)
		sr1, dr1 := vadd(r1, r3), vsub(r1, r3)
		si0, di0 := vadd(i0, i2), vsub(i0, i2)
		si1, di1 := vadd(i1, i3), vsub(i1, i3)

		r0, i0 = vadd(sr0, sr1), vadd(si0, si1)
		r1, i1 = vsub(dr0, di1), vadd(di0, dr1)
		r2, i2 = vsub(sr0, sr1), vsub(si0, si1)
		r3, i3 = vadd(dr0, di1), vsub(di0, dr1)

		r1, i1 = vcplxmulconj(r1, i1, e[k*6+0], e[k*6+1])
		r2, i2 = vcplxmulconj(r2, i2, e[k*6+2], e[k*6+3])
		r3, i3 = vcplxmulconj(r3, i3, e[k*6+4], e[k*6+5])

		r0, r1, r2, r3 = vtranspose4(r0, r1, r2, r3)
		i0, i1, i2, i3 = vtranspose4(i0, i1, i2, i3)

		out[0], out[1], out[2], out[3] = r0, i0, r1, i1
		out[4], out[5], out[6], out[7] = r2, i2, r3, i3
		out = out[8:]
	}
}

func pffft_real_finalize_4x4(in0, in1 reg, in, e, out []reg) {
	r0, i0 := in0, in1
	r1, i1, r2, i2, r3, i3 := in[0], in[1], in[2], in[3], in[4], in[5]
	r0, r1, r2, r3 = vtranspose4(r0, r1, r2, r3)
	i0, i1, i2, i3 = vtranspose4(i0, i1, i2, i3)

	/*
	   transformation for each column is:

	   [1   1   1   1   0   0   0   0]   [r0]
	   [1   0  -1   0   0  -1   0   1]   [r1]
	   [1   0  -1   0   0   1   0  -1]   [r2]
	   [1  -1   1  -1   0   0   0   0]   [r3]
	   [0   0   0   0   1   1   1   1] * [i0]
	   [0  -1   0   1  -1   0   1   0]   [i1]
	   [0  -1   0   1   1   0  -1   0]   [i2]
	   [0   0   0   0  -1   1  -1   1]   [i3]
	*/
	r1, i1 = vcplxmul(r1, i1, e[0], e[1])
	r2, i2 = vcplxmul(r2, i2, e[2], e[3])
	r3, i3 = vcplxmul(r3, i3, e[4], e[5])

	sr0, dr0 := vadd(r0, r2), vsub(r0, r2)
	sr1, dr1 := vadd(r1, r3), vsub(r3, r1)
	si0, di0 := vadd(i0, i2), vsub(i0, i2)
	si1, di1 := vadd(i1, i3), vsub(i3, i1)

	r0 = vadd(sr0, sr1)
	r3 = vsub(sr0, sr1)
	i0 = vadd(si0, si1)
	i3 = vsub(si1, si0)
	r1 = vadd(dr0, di1)
	r2 = vsub(dr0, di1)
	i1 = vsub(dr1, di0)
	i2 = vadd(dr1, di0)

	out[0] = r0
	out[1] = i0
	out[2] = r1
	out[3] = i1
	out[4] = r2
	out[5] = i2
	out[6] = r3
	out[7] = i3
}

//go:noinline
func pffft_real_finalize(ncvec int, in, out, e []reg) {
	SIMD_SZ := simd.VectorBitSize() / 32
	dk := ncvec / SIMD_SZ // number of 4x4 matrix blocks
	if &in[0] == &out[0] {
		panic(`gopffft: pffft_cplx_preprocess: in and out can't alias`)
	}

	var cr, ci [32]float32
	save := in[0]
	const s = float32(1 / math.Sqrt2)

	in[0].Store(cr[:])
	in[ncvec*2-1].Store(ci[:])
	pffft_real_finalize_4x4(reg{}, reg{}, in[1:], e, out)

	/*
	   [cr0 cr1 cr2 cr3 ci0 ci1 ci2 ci3]

	   [Xr(1)]  ] [1   1   1   1   0   0   0   0]
	   [Xr(N/4) ] [0   0   0   0   1   s   0  -s]
	   [Xr(N/2) ] [1   0  -1   0   0   0   0   0]
	   [Xr(3N/4)] [0   0   0   0   1  -s   0   s]
	   [Xi(1)   ] [1  -1   1  -1   0   0   0   0]
	   [Xi(N/4) ] [0   0   0   0   0  -s  -1  -s]
	   [Xi(N/2) ] [0  -1   0   1   0   0   0   0]
	   [Xi(3N/4)] [0   0   0   0   0  -s   1  -s]
	*/

	out[0].StorePart([]float32{(cr[0] + cr[2]) + (cr[1] + cr[3])}) // xr0
	out[1].StorePart([]float32{(cr[0] + cr[2]) - (cr[1] + cr[3])}) // xi0
	out[4].StorePart([]float32{(cr[0] - cr[2])})                   // xr2
	out[5].StorePart([]float32{(cr[3] - cr[1])})                   // xi2
	out[2].StorePart([]float32{ci[0] + s*(ci[1]-ci[3])})           // xr1
	out[3].StorePart([]float32{-ci[2] - s*(ci[1]+ci[3])})          // xi1
	out[6].StorePart([]float32{ci[0] - s*(ci[1]-ci[3])})           // xr3
	out[7].StorePart([]float32{ci[2] - s*(ci[1]+ci[3])})           // xi3

	for k := 1; k < dk; k++ {
		save_next := in[8*k+7]
		pffft_real_finalize_4x4(save, in[8*k+0], in[8*k+1:], e[k*6:], out[k*8:])
		save = save_next
	}
}

func pffft_real_preprocess_4x4(in, e, out []reg, first bool) {
	r0, i0, r1, i1, r2, i2, r3, i3 := in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7]
	/*
	   transformation for each column is:

	   [1   1   1   1   0   0   0   0]   [r0]
	   [1   0   0  -1   0  -1  -1   0]   [r1]
	   [1  -1  -1   1   0   0   0   0]   [r2]
	   [1   0   0  -1   0   1   1   0]   [r3]
	   [0   0   0   0   1  -1   1  -1] * [i0]
	   [0  -1   1   0   1   0   0   1]   [i1]
	   [0   0   0   0   1   1  -1  -1]   [i2]
	   [0   1  -1   0   1   0   0   1]   [i3]
	*/
	sr0, dr0 := vadd(r0, r3), vsub(r0, r3)
	sr1, dr1 := vadd(r1, r2), vsub(r1, r2)
	si0, di0 := vadd(i0, i3), vsub(i0, i3)
	si1, di1 := vadd(i1, i2), vsub(i1, i2)

	r0 = vadd(sr0, sr1)
	r2 = vsub(sr0, sr1)
	r1 = vsub(dr0, si1)
	r3 = vadd(dr0, si1)
	i0 = vsub(di0, di1)
	i2 = vadd(di0, di1)
	i1 = vsub(si0, dr1)
	i3 = vadd(si0, dr1)

	r1, i1 = vcplxmulconj(r1, i1, e[0], e[1])
	r2, i2 = vcplxmulconj(r2, i2, e[2], e[3])
	r3, i3 = vcplxmulconj(r3, i3, e[4], e[5])

	r0, r1, r2, r3 = vtranspose4(r0, r1, r2, r3)
	i0, i1, i2, i3 = vtranspose4(i0, i1, i2, i3)

	if !first {
		out[0] = r0
		out[1] = i0
		out = out[2:]
	}
	out[0] = r1
	out[1] = i1
	out[2] = r2
	out[3] = i2
	out[4] = r3
	out[5] = i3
}

//go:noinline
func pffft_real_preprocess(ncvec int, in, out, e []reg) {
	SIMD_SZ := simd.VectorBitSize() / 32
	dk := ncvec / SIMD_SZ // number of 4x4 matrix blocks
	if &in[0] == &out[0] {
		panic(`gopffft: pffft_cplx_preprocess: in and out can't alias`)
	}

	var Xr, Xi [32]float32
	const s = float32(math.Sqrt2)
	for k := range 4 {
		var t [32]float32
		in[(8*k)/SIMD_SZ].Store(t[:])
		Xr[k] = t[0]
		in[(8*k+4)/SIMD_SZ].Store(t[:])
		Xi[k] = t[0]
	}

	pffft_real_preprocess_4x4(in, e, out[1:], true) // will write only 6 values
	/*
	   [Xr0 Xr1 Xr2 Xr3 Xi0 Xi1 Xi2 Xi3]

	   [cr0] [1   0   2   0   1   0   0   0]
	   [cr1] [1   0   0   0  -1   0  -2   0]
	   [cr2] [1   0  -2   0   1   0   0   0]
	   [cr3] [1   0   0   0  -1   0   2   0]
	   [ci0] [0   2   0   2   0   0   0   0]
	   [ci1] [0   s   0  -s   0  -s   0  -s]
	   [ci2] [0   0   0   0   0  -2   0   2]
	   [ci3] [0  -s   0   s   0  -s   0  -s]
	*/
	for k := 1; k < dk; k++ {
		pffft_real_preprocess_4x4(in[8*k:], e[k*6:], out[-1+k*8:], false)
	}
	cr0 := (Xr[0] + Xi[0]) + 2*Xr[2]
	cr1 := (Xr[0] - Xi[0]) - 2*Xi[2]
	cr2 := (Xr[0] + Xi[0]) - 2*Xr[2]
	cr3 := (Xr[0] - Xi[0]) + 2*Xi[2]
	out[0].StorePart([]float32{cr0, cr1, cr2, cr3})
	ci0 := 2 * (Xr[1] + Xr[3])
	ci1 := s*(Xr[1]-Xr[3]) - s*(Xi[1]+Xi[3])
	ci2 := 2 * (Xi[3] - Xi[1])
	ci3 := -s*(Xr[1]-Xr[3]) - s*(Xi[1]+Xi[3])
	out[2*ncvec-1].StorePart([]float32{ci0, ci1, ci2, ci3})
}

func pffft_transform_internal(s *pffftSetup, finput, foutput []float32, scratch []reg, forward bool, ordered int) (ainput, aoutput []float32) {
	ncvec := s.ncvec
	nf_odd := s.ifac[1] & 1

	// NOTE(neputevshina): No stack allocation.

	ainput = align(finput)
	aoutput = align(foutput)

	vinput := getAsAligned(ainput)
	voutput := getAsAligned(aoutput)
	buff := [2][]reg{voutput, scratch}

	ib := nf_odd ^ ordered

	if forward {
		ib = 1 - ib
	} else {
		if &vinput[0] == &buff[ib][0] {
			ib = 1 - ib // may happen when finput == foutput
		}
		if ordered == 1 {
			pffft_zreorder(s, degrade(vinput), degrade(buff[ib]), false)
			vinput = buff[ib]
			ib = 1 - ib
		}
	}
	switch {
	case forward && s.transform == transformReal:
		o := rfftf1_ps(ncvec*2, vinput, buff[ib], buff[1-ib], s.twiddle, s.ifac[:])
		if &o[0] == &buff[0][0] {
			ib = 0
		} else {
			ib = 1
		}
		pffft_real_finalize(ncvec, buff[ib], buff[1-ib], s.e)
	case forward && s.transform == transformComplex:
		tmp := buff[ib]
		for k := range ncvec {
			tmp[k*2], tmp[k*2+1] = uninterleave2(vinput[k*2], vinput[k*2+1])
		}
		o := cfftf1_ps(ncvec, buff[ib], buff[1-ib], buff[ib], s.twiddle, s.ifac[:], -1)
		if &o[0] == &buff[0][0] {
			ib = 0
		} else {
			ib = 1
		}
		pffft_cplx_finalize(ncvec, buff[ib], buff[1-ib], s.e)
	case !forward && s.transform == transformReal:
		pffft_real_preprocess(ncvec, vinput, buff[ib], s.e)
		o := rfftb1_ps(ncvec*2, buff[ib], buff[0], buff[1], s.twiddle, s.ifac[:])
		if &o[0] == &buff[0][0] {
			ib = 0
		} else {
			ib = 1
		}
	case !forward && s.transform == transformComplex:
		pffft_cplx_preprocess(ncvec, vinput, buff[ib], s.e)
		o := cfftf1_ps(ncvec, buff[ib], buff[0], buff[1], s.twiddle, s.ifac[:], +1)
		if &o[0] == &buff[0][0] {
			ib = 0
		} else {
			ib = 1
		}
		for k := range ncvec {
			buff[ib][k*2], buff[ib][k*2+1] = interleave2(buff[ib][k*2], buff[ib][k*2+1])
		}
	}

	if &buff[ib][0] != &voutput[0] {
		/* extra copy required -- this situation should only happen when finput == foutput */
		if &finput[0] != &foutput[0] {
			panic(`incorrect`)
		}
		for k := range ncvec {
			a, b := buff[ib][2*k], buff[ib][2*k+1]
			voutput[2*k], voutput[2*k+1] = a, b
		}
		ib = 1 - ib
	}
	if &buff[ib][0] != &voutput[0] {
		panic(`incorrect`)
	}

	return
}

func pffft_transform(s *pffftSetup, input, output, work []float32, forward bool) (ainput, aoutput []float32) {
	return pffft_transform_internal(s, input, output, getAsAligned(work), forward, 0)
}

func pffft_transform_ordered(s *pffftSetup, input, output, work []float32, forward bool) (ainput, aoutput []float32) {
	return pffft_transform_internal(s, input, output, getAsAligned(work), forward, 1)
}
