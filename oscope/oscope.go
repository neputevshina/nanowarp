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
	"strconv"

	"golang.org/x/exp/constraints"
)

var Enable bool

type where [2]uintptr

var watches = map[where]*watch{}

type Etc any

type watch struct {
	outfile string
	elem    any
	data    []any
	etc     []Etc
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
		watches[wh] = &watch{
			outfile: fmt.Sprintf("%s:%d(%d)", path.Base(file), line, wh[1]),
			etc:     etc,
			elem:    a,
		}
		w = watches[wh]
	}
	w.data = append(w.data, a)
}

func Dump(err error, topath string) error {
	if err != nil {
		return err
	}
	if !Enable {
		return nil
	}
	for _, w := range watches {
		outf := path.Join(topath, w.outfile+".png")
		switch w.elem.(type) {
		case float64:
			err = dumpWaveform[float64](nil, w, w.data, outf)
		case float32:
			err = dumpWaveform[float32](nil, w, w.data, outf)
		case int:
			err = dumpWaveform[int](nil, w, w.data, outf)
		case uint8:
			err = dumpWaveform[uint8](nil, w, w.data, outf)

		case []float64:
			err = dumpTexture[float64](nil, w, w.data, outf)
		case []float32:
			err = dumpTexture[float32](nil, w, w.data, outf)
		case []int:
			err = dumpTexture[int](nil, w, w.data, outf)
		case []uint8:
			err = dumpTexture[uint8](nil, w, w.data, outf)
		}
		if err != nil {
			return err
		}
	}
	return nil
}

func dumpWaveform[T constraints.Integer | constraints.Float](err error, w *watch, data []any, topath string) error {
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
	for x := 0; x < width; x++ {
		y := float64(data[x].(T))*scale + offset
		img.SetGray(x, int(math.Floor(y)), color.Gray{Y: 255})
	}

	file, err := os.Create(topath)
	if err != nil {
		return err
	}
	return png.Encode(file, img)
}

func dumpTexture[T constraints.Integer | constraints.Float](err error, w *watch, data []any, topath string) error {
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
	for x := 0; x < width; x++ {
		for y := 0; y < height; y++ {
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
