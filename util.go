package nanowarp

import (
	"fmt"
	"image"
	"image/color"
	"math"
	"math/cmplx"
	"reflect"

	"golang.org/x/exp/constraints"
)

var mag = cmplx.Abs

type bang = struct{}

var println = fmt.Println

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

func clamp[T constraints.Ordered](a, b, x T) T {
	return max(a, min(b, x))
}

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

func getfadv(x, xt []complex128, stretch float64) func(w int) float64 {
	return func(j int) float64 {
		if mag(x[j]) < 1e-6 {
			return 0
		}
		// TODO This phase correction value is guaranteed to be wrong but is mostly correct.
		return -real(xt[j]/x[j])/float64(len(x))*math.Pi*stretch - math.Pi/2
	}
}

type heaptriple struct {
	mag  float64
	w, t int
}
type hp []heaptriple

func (h hp) Len() int           { return len(h) }
func (h hp) Less(i, j int) bool { return h[i].mag > h[j].mag }
func (h hp) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *hp) Push(x any) {
	*h = append(*h, x.(heaptriple))
}
func (h *hp) Pop() any {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func hann(out []float64) {
	for i := range out {
		x := float64(i) / float64(len(out))
		out[i] = 0.5 * (1.0 - math.Cos(2.0*math.Pi*x))
	}
}

func hannDx(out []float64) {
	for i := range out {
		x := float64(i)/float64(len(out)) + .5
		out[i] = math.Pi * math.Sin(2*math.Pi*x)
	}
}

func niemitalo(out []float64) {
	// https://dsp.stackexchange.com/questions/2337/fft-with-asymmetric-windowing
	nfft := float64(len(out))
	clear(out)
	sin, cos := math.Sin, math.Cos
	for i := nfft / 4; i < nfft*7/8; i++ {
		x := 2 * math.Pi * ((i+0.5)/nfft - 1.75)
		out[int(i)] = 2.57392230162633461887 - 1.58661480271141974718*cos(x) + 3.80257516644523141380*sin(x) -
			1.93437090055110760822*cos(2*x) - 3.27163999159752183488*sin(2*x) + 3.26617449847621266201*cos(3*x) -
			0.30335261753524439543*sin(3*x) - 0.92126091064427817479*cos(4*x) + 2.33100177294084742741*sin(4*x) -
			1.19953922321306438725*cos(5*x) - 1.25098147932225423062*sin(5*x) + 0.99132076607048635886*cos(6*x) -
			0.34506787787355830410*sin(6*x) - 0.04028033685700077582*cos(7*x) + 0.55461815542612269425*sin(7*x) -
			0.21882110175036428856*cos(8*x) - 0.10756484378756643594*sin(8*x) + 0.06025986430527170007*cos(9*x) -
			0.05777077835678736534*sin(9*x) + 0.00920984524892982936*cos(10*x) + 0.01501989089735343216*sin(10*x)
	}
	for i := 0; i < int(nfft)/8; i++ {
		nfft := int(nfft)
		out[nfft-1-i] = (1 - out[nfft*3/4-1-i]*out[nfft*3/4+i]) / out[nfft/2+i]
	}
	copy(out, out[int(nfft)*2/8:])
	clear(out[int(nfft)*6/8:])
}

func niemitaloDx(out []float64) {
	panic(`unimplemented`)
}

func windowGain(w []float64) (a float64) {
	for _, e := range w {
		a += e * e
	}
	a /= float64(len(w))
	return
}

func windowT(w, out []float64) {
	n := float64(len(w))
	for i := range w {
		out[i] = w[i] * mix(-n/2, n/2+1, float64(i)/n)
	}
}

func floatMatrixToImage(data [][]float64) image.Image {
	if len(data) == 0 || len(data[0]) == 0 {
		return nil
	}

	height := len(data[0])
	width := len(data)

	// Find min and max
	minVal := math.Inf(1)
	maxVal := math.Inf(-1)
	for _, row := range data {
		for _, v := range row {
			if v < minVal {
				minVal = v
			}
			if v > maxVal {
				maxVal = v
			}
		}
	}
	fmt.Println(minVal, maxVal)

	scale := 1.
	offset := 3.14
	scale = 255.0 / (maxVal - minVal)
	offset = -minVal * scale

	img := image.NewGray(image.Rect(0, 0, width, height))
	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			val := data[x][y]*scale + offset
			if val < 0 {
				val = 0
			}
			if val > 255 {
				val = 255
			}
			img.SetGray(x, y, color.Gray{Y: uint8(val + 0.5)})
		}
	}

	return img
}
