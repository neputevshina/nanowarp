package oscope

import (
	"bytes"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
	"path"
	"runtime"
	"slices"
	"strconv"

	"golang.org/x/exp/constraints"
)

var Enable bool

type where [2]uintptr

var watches = map[where]*watch{}

type Etc = any

type Name string
type Normalize bool
type Quantization int

type watch struct {
	outfile string
	elem    any
	data    []any
	etc     []Etc
	norm    bool
}

// It's okay to get goroutine ID for debugging purposes,
// since this is a graphical debugging package.
func getGID() uintptr {
	b := make([]byte, 64)
	b = b[:runtime.Stack(b, false)]
	b = bytes.TrimPrefix(b, []byte("goroutine "))
	b = b[:bytes.IndexByte(b, ' ')]
	n, _ := strconv.ParseUint(string(b), 10, 64)
	return uintptr(n)
}

func findetc[T any](etc []Etc) (int, T, bool) {
	i := slices.IndexFunc(etc, func(e Etc) bool { _, k := e.(T); return k })
	var z T
	if i >= 0 {
		z = etc[i].(T)
	}
	return i, z, i >= 0
}

func Oscope(a any, etc ...Etc) {
	if !Enable {
		return
	}
	pc, file, line, ok := runtime.Caller(1)
	wh := [2]uintptr{pc, getGID()}
	w, k := watches[wh]
	if !k {
		if !ok {
			panic(`oscope.Oscope: can't identify a function that is not traceable on the stack`)
		}
		fn := ""
		_, e, ok := findetc[Name](etc)
		if ok {
			fn = string(e)
		} else {
			fn = fmt.Sprintf("%s:%d(%d)", path.Base(file), line, wh[1])
		}

		watches[wh] = &watch{
			outfile: fn,
			etc:     etc,
			elem:    a,
		}
		w = watches[wh]
		_, n, _ := findetc[Normalize](etc)
		w.norm = bool(n)
	}
	w.data = append(w.data, a)
}

func dumpWaveform[T constraints.Integer | constraints.Float](err error, w *watch, topath string) error {
	topath = path.Join(topath, w.outfile+".png")
	data := w.data
	n, x := data[0].(T), data[0].(T)
	for _, e := range data {
		n = min(n, e.(T))
		x = max(x, e.(T))
	}

	height := 1024.
	width := len(data)

	scale := height / float64(x-n)
	offset := float64(-n) * scale

	img := image.NewGray(image.Rect(0, 0, width, int(height)))
	for x := range width {
		y := 1 - (float64(data[x].(T))*scale + offset)
		img.SetGray(x, int(math.Floor(y)), color.Gray{Y: 255})
	}

	file, err := os.Create(topath)
	if err != nil {
		return err
	}
	return png.Encode(file, img)
}

func dumpTexture[T constraints.Integer | constraints.Float](err error, w *watch, topath string) error {
	topath = path.Join(topath, w.outfile+".png")
	data := w.data
	type S = []T
	n, x := data[0].(S)[0], data[0].(S)[0]
	for _, s := range data {
		for _, e := range s.(S) {
			n = min(n, e)
			x = max(x, e)
		}
	}

	height := len(data[0].(S))
	width := len(data)

	scale := 255.0 / float64(x-n)
	offset := float64(-n) * scale

	img := image.NewGray(image.Rect(0, 0, width, height))
	for x := range width {
		for y := range height {
			v := float64(data[x].(S)[y])*scale + offset
			img.SetGray(x, y, color.Gray{Y: uint8(max(0, min(255, v+0.5)))})
		}
	}

	file, err := os.Create(topath)
	if err != nil {
		return err
	}
	return png.Encode(file, img)
}
