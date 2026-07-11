package nanowarp

import (
	"cmp"
	"fmt"
	"math"
	"math/cmplx"
	"reflect"
	"strconv"

	"golang.org/x/exp/constraints"
)

var println = fmt.Println

type bang = struct{}

func bitsafe(v float64) float64 {
	if v != v || math.IsInf(v, 0) {
		return 0
	}
	return v
}

func princarg(phase float64) float64 {
	pi2 := 2 * math.Pi
	return phase - math.Round(phase/pi2)*pi2
}

func add[T constraints.Float | constraints.Complex](dst, src []T) []T {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] += src[i]
	}
	return dst
}

func sum[T constraints.Float | constraints.Complex](src []T) (s T) {
	for i := range src {
		s += src[i]
	}
	return
}

func mul[T constraints.Float | constraints.Complex](dst, src []T) []T {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] *= src[i]
	}
	return dst
}

// mix is a linear interpolation.
func mix[F constraints.Float](a, b, x F) F {
	return a*(1-x) + b*x
}

// precisionmix is a highly accurate linear interpolation.
// It is slow because math.FMA is done in software floating point.
func precisionmix(a, b, x float64) float64 {
	return math.FMA(x, b-a, a)
}

// unmix is a linear extrapolation.
func unmix[F constraints.Float](a, b, x F) F {
	return (x - a) / (b - a)
}

// clamp hard clips value x to satisfy a ≤ x ≤ b.
func clamp[T cmp.Ordered](a, b, x T) T {
	return max(a, min(b, x))
}

// makeslices is a helper function to initialize all slices in a
// structure instance using reflection.
//
// a must be a pointer to struct.
func makeslices(a any, nbins, nfft, nch, lah int) {
	rn := reflect.ValueOf(a).Elem()
	tf := reflect.TypeOf(a).Elem()
	for i := 0; i < rn.NumField(); i++ {
		f := rn.Field(i)
		ln := map[string]int{
			`nbins`: nbins,
			`nfft`:  nfft,
			`lah`:   lah,
			`nch`:   nch,
		}[tf.Field(i).Tag.Get("size")]

		if f.Kind() == reflect.Slice {
			c := f.Interface()
			switch c.(type) {
			case []complex128:
				if ln == 0 {
					ln = nbins
				}
				f.Set(reflect.ValueOf(make([]complex128, ln)))
			case []float64:
				if ln == 0 {
					ln = nfft
				}
				f.Set(reflect.ValueOf(make([]float64, ln)))
			case [][]complex128:
				if ln == 0 {
					ln = nch
				}
				f.Set(reflect.ValueOf(make([][]complex128, ln)))
				s := f.Interface().([][]complex128)
				for i := range s {
					s[i] = make([]complex128, nbins)
				}
			case [][][]complex128:
				f.Set(reflect.ValueOf(make([][][]complex128, lah)))
				s := f.Interface().([][][]complex128)
				for ch := range s {
					s[ch] = make([][]complex128, nch)
					for i := range s[ch] {
						s[ch][i] = make([]complex128, nbins)
					}
				}
			case [][]float64:
				f.Set(reflect.ValueOf(make([][]float64, ln)))
				s := f.Interface().([][]float64)
				for i := range s {
					s[i] = make([]float64, nbins)
				}
			}
		}
	}
}

// structinit replaces all zero values in a given struct with a value
// from `default` annotation for each corresponding field.
//
// a must be a pointer to struct.
func structinit(a any) {
	rn := reflect.ValueOf(a).Elem()
	ts := reflect.TypeOf(a).Elem()
	for i := 0; i < rn.NumField(); i++ {
		f := rn.Field(i)
		if !f.IsZero() {
			continue
		}
		tf := ts.Field(i)
		v := tf.Tag.Get("default")
		if v == `` {
			continue
		}
		switch tf.Type {
		case reflect.TypeFor[int]():
			vv, err := strconv.ParseInt(v, 10, 64)
			if err != nil {
				panic(fmt.Sprint(`incorrect default value on field `, tf.Name, ` in type `, ts.Name()))
			}
			f.Set(reflect.ValueOf(int(vv)))
		case reflect.TypeFor[float64]():
			vv, err := strconv.ParseFloat(v, 64)
			if err != nil {
				panic(fmt.Sprint(`incorrect default value on field `, tf.Name, `in type `, ts.Name()))
			}
			f.Set(reflect.ValueOf(vv))
		case reflect.TypeFor[string]():
			f.Set(reflect.ValueOf(v))
		}

	}
}

// nextpow2 returns the minimum power of two greater than i.
func nextpow2(i int) int {
	return int(math.Floor(math.Pow(2, math.Ceil(math.Log2(float64(i))))))
}

func norm(c complex128) complex128 {
	m := cmplx.Abs(c)
	if m == 0 {
		return 0
	}
	return c / complex(cmplx.Abs(c), 0)
}

func fill[T any](s []T, e T) {
	for i := range s {
		s[i] = e
	}
}

func abs[T constraints.Signed | constraints.Float](a T) T {
	if a < 0 {
		return -a
	}
	return a
}

func boolfloat(b bool) float64 {
	if b {
		return 1
	}
	return 0
}

// hztobin converts a frequency value in hertz to bin number of a DFT
// of size nfft on a signal with known sample rate.
func hztobin(hz float64, nfft, samplerate int) int {
	return int(hz * float64(nfft) / float64(samplerate))
}

func scale[T constraints.Float | constraints.Complex](dst []T, s T) {
	for i := range dst {
		dst[i] *= s
	}
}

func even(x int) int {
	return x + x%2
}

func softmax(a []float64) {
	expsum := 0.
	for _, v := range a {
		expsum += math.Exp2(v)
	}
	for i := range a {
		a[i] = math.Exp2(a[i]) / expsum
	}
}

func make2[T any](j, i int) (v [][]T) {
	v = make([][]T, j)
	for j := range j {
		v[j] = make([]T, i)
	}
	return
}

func make3[T any](k, j, i int) (v [][][]T) {
	v = make([][][]T, j)
	for k := range k {
		v[k] = make([][]T, i)
		for j := range j {
			v[k][j] = make([]T, i)
		}
	}
	return
}

func safediv[T constraints.Float | constraints.Complex](a, b T) T {
	if b == 0 {
		return 0
	}
	return a / b
}

// atodb converts an amplidude value to decibels full scale (dBFS).
func atodb(a float64) float64 {
	return 20 * math.Log10(a)
}

// dbtoa converts a value in decibels full scale (dBFS) to amplitude.
func dbtoa(db float64) float64 {
	return math.Pow(10, db/20)
}
