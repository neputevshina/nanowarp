package dspio

import (
	"errors"
	"fmt"
	"io"
	"sync/atomic"

	"golang.org/x/sys/cpu"
)

var ErrSpillover error = errors.New(`buffer spillover: wait, lock or timeout`)

// Inspired by
// https://blog.systems.ethz.ch/blog/2019/the-design-and-implementation-of-a-lock-free-ring-buffer-with-contiguous-reservations.html

type altpipe struct {
	bufs [][]float64
	nch  int

	_     cpu.CacheLinePad
	read  atomic.Int64
	_     cpu.CacheLinePad
	write atomic.Int64
	_     cpu.CacheLinePad
}

func (p *altpipe) NchWrite() int { return p.nch }
func (p *altpipe) NchRead() int  { return p.nch }

func (p *altpipe) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(buf) != p.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, p.nch, len(buf)))
	}

	if int(p.write.Load()) == len(p.bufs[0]) {
		p.write.Store(0)
	}
	if p.write.Load() == p.read.Load() {
		// Queue is full.
		return 0, ErrSpillover
	}
	r, w := int(p.read.Load()), int(p.write.Load())
	if w < r {
		n = copy(p.bufs[0][w:r], buf[0])
		p.write.Add(int64(n))
	} else {
		a := copy(p.bufs[0][w:], buf[0])
		b := copy(p.bufs[0][:r], buf[0][a:])
		n = a + b
		p.write.Store(int64(b))
	}
	if n < len(buf[0]) {
		err = io.ErrShortWrite
	}

	return
}

func (p *altpipe) SignalRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(buf) != p.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, p.nch, len(buf)))
	}

	if int(p.read.Load()) == len(p.bufs[0]) {
		p.read.Store(0)
	}
	// NOTE We can always read, if p.read == p.write we'll just copy nothing and n will be 0.
	r, w := int(p.read.Load()), int(p.write.Load())
	if w < r {
		a := copy(buf[0], p.bufs[0][r:])
		b := copy(buf[0][a:], p.bufs[0][:w])
		n = a + b
		p.read.Store(int64(b))
	} else {
		n = copy(buf[0], p.bufs[0][r:w])
		p.read.Add(int64(n))
	}

	return
}

var _ SignalReader = &altpipe{}
var _ SignalWriter = &altpipe{}

func Pipe(nch int, nbuf int) (SignalWriter, SignalReader) {
	p := &altpipe{
		nch:  nch,
		bufs: make2(nch, nbuf),
	}

	return p, p
}

func Copy(prr error, r SignalReader, w SignalWriter, buf [][]float64) (written int, err error) {
	if prr != nil {
		return 0, prr
	}
	if buf == nil {
		buf = make2(r.NchRead(), 8192)
	}
	if len(buf) != r.NchRead() && len(buf) != w.NchWrite() {
		panic(fmt.Errorf(`different number of channels: got r=%d, buf=%d, w=%d`, r.NchRead(), w.NchWrite(), len(buf)))
	}

	return
}
