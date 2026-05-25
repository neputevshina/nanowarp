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
//
// DOES NOT WORK

type lockfreePipeWriter struct {
	bufs [][]float64
	nch  int
	_    cpu.CacheLinePad

	read  atomic.Int64
	_     cpu.CacheLinePad
	write atomic.Int64
	_     cpu.CacheLinePad
	close atomic.Bool
	_     cpu.CacheLinePad
}

type lockfreePipeReader struct {
	*lockfreePipeWriter
}

func (p *lockfreePipeWriter) NchWrite() int { return p.nch }
func (p *lockfreePipeReader) NchRead() int  { return p.nch }

func (p *lockfreePipeWriter) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if p.close.Load() {
		return 0, io.EOF
	}
	if len(buf) != p.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, p.nch, len(buf)))
	}

	if int(p.write.Load()) == len(p.bufs[0]) {
		p.write.Store(0)
	}
	if p.write.Load() == p.read.Load()-1 {
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
		if a+w < len(p.bufs[0]) {
			p.write.Store(int64(n))
		} else {
			p.write.Store(int64(b))
		}
	}
	if n < len(buf[0]) {
		err = io.ErrShortWrite
	}

	return
}

func (p *lockfreePipeReader) SignalRead(prr error, buf [][]float64) (n int, err error) {
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
	if p.close.Load() && n == 0 {
		return 0, io.EOF
	}

	return
}

func (p *lockfreePipeWriter) Close() (err error) {
	p.close.Store(true)
	return nil
}

var _ SignalWriter = &lockfreePipeWriter{}
var _ SignalReader = &lockfreePipeReader{}

func lockfreePipe(nch int, nbuf int) (*lockfreePipeWriter, *lockfreePipeReader) {
	p := &lockfreePipeWriter{
		nch:  nch,
		bufs: make2(nch, nbuf),
	}

	return p, &lockfreePipeReader{p}
}

func Copy(prr error, r SignalReader, w SignalWriter, buf [][]float64) (written int, err error) {
	if prr != nil {
		return 0, prr
	}
	if buf == nil {
		buf = make2(r.NchRead(), 8192)
	}
	sl := make([][]float64, r.NchRead())

	n, no := 0, 0
	defer func() {
		if err == io.EOF {
			err = nil
		}
	}()
	for err == nil {
		n, err = r.SignalRead(nil, buf)
		if n <= 0 {
			continue
		}
		if err != nil {
			return
		}
		for i := range sl {
			sl[i] = buf[i][:n]
		}
		no, err = w.SignalWrite(nil, sl)
		if no < 0 || n < no {
			if err == nil {
				err = fmt.Errorf(`dspio.Copy: invalid write result, written %v, expected to write %v`, no, n)
			}
			return
		}
		written += no
	}

	return
}
