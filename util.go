package nanowarp

import (
	"math"
	"math/cmplx"
	"reflect"

	"golang.org/x/exp/constraints"
)

// G is an implicit global variable map for internal debugging purposes.
// Accesses to this map are either in a thing that is not done yet or safe to remove.
//
// You MUST NOT initialize this map.
// If something fails because of it, it's because Nanowarp is broken,
// and you must use a previous version of the library.
var G map[string]any

var mag = cmplx.Abs

type bang = struct{}

func princarg(phase float64) float64 {
	pi2 := 2 * math.Pi
	return phase - math.Round(phase/pi2)*pi2
}

func add[T constraints.Float](dst, src []T) {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] += src[i]
	}
}

func sub[T constraints.Float](dst, src []T) {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] -= src[i]
	}
}

func mul[T constraints.Float](dst, src []T) {
	for i := 0; i < min(len(dst), len(src)); i++ {
		dst[i] *= src[i]
	}
}

func mix[F constraints.Float](a, b, x F) F {
	return a*(1-x) + b*x
}

// makeslices automatically initializes slices through reflection.
func makeslices(a any, nbins, nfft int) {
	rn := reflect.ValueOf(a).Elem()
	for i := 0; i < rn.NumField(); i++ {
		f := rn.Field(i)
		if f.Kind() == reflect.Slice {
			c := f.Interface()
			switch c := c.(type) {
			case []complex128:
				f.Set(reflect.ValueOf(make([]complex128, nbins)))
			case []float64:
				f.Set(reflect.ValueOf(make([]float64, nfft)))
			case [][]complex128:
				for i := range c {
					c[i] = make([]complex128, nbins)
				}
			}
		}
	}
}
