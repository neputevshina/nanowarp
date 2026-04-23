package dspio

import (
	"errors"
	"io"
)

type SignalReader interface {
	NchRead() int
	SignalRead(prr error, buf [][]float64) (n int, err error)
}

type SignalWriter interface {
	NchWrite() int
	SignalWrite(prr error, buf [][]float64) (n int, err error)
}

func ReadAll(prr error, r SignalReader) (d [][]float64, err error) {
	if prr != nil {
		return nil, prr
	}

	nch := r.NchRead()
	t := make2(nch, 8192)
	d = make([][]float64, nch)
	n := 0

	for err == nil {
		n, err = r.SignalRead(nil, t)
		for ch := range nch {
			d[ch] = append(d[ch], t[ch][:n]...)
		}
	}

	if err == io.EOF {
		err = nil
	}

	return
}

var Discard SignalWriter = &discard{}

type discard struct{}

func (d *discard) NchWrite() int { return 0 }

func (d *discard) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	if prr == nil {
		n = len(buf[0])
	}
	return n, prr
}

func TeeReader(r SignalReader, w SignalWriter) SignalReader {
	return &teereader{r, w, make([][]float64, r.NchRead())}
}

type teereader struct {
	r     SignalReader
	w     SignalWriter
	knife [][]float64
}

func (t *teereader) NchRead() int { return t.r.NchRead() }

func (t *teereader) SignalRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	n, err = t.r.SignalRead(nil, buf)
	if n > 0 {
		for ch := range t.r.NchRead() {
			t.knife[ch] = buf[ch][:n]
		}
		if n, err := t.w.SignalWrite(nil, t.knife); err != nil {
			return n, err
		}
	}
	return
}

func LimitReader(r SignalReader, sa int) SignalReader {
	return &limitreader{sa, r, make([][]float64, r.NchRead())}
}

type limitreader struct {
	N     int
	R     SignalReader
	knife [][]float64
}

func (t *limitreader) NchRead() int { return t.R.NchRead() }

func (l *limitreader) SignalRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if l.N <= 0 {
		return 0, io.EOF
	}
	for ch := range l.R.NchRead() {
		if len(buf[0]) > l.N {
			l.knife[ch] = buf[ch][:l.N]
		} else {
			l.knife[ch] = buf[ch]
		}
	}
	n, err = l.R.SignalRead(nil, l.knife)
	l.N -= n
	return
}

var ErrNonSerialWrite error = errors.New(`non-serial write: amount of samples written is not equal for all channels`)
var ErrNonSerialRead error = errors.New(`non-serial read: amount of samples consumed is not equal for all channels`)
