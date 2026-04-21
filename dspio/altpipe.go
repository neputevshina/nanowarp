package dspio

import (
	"fmt"
	"sync/atomic"

	"golang.org/x/sys/cpu"
)

// https://blog.systems.ethz.ch/blog/2019/the-design-and-implementation-of-a-lock-free-ring-buffer-with-contiguous-reservations.html
type altpipe struct {
	bufs [][]float64
	nch  int

	_         cpu.CacheLinePad
	read      atomic.Int64
	_         cpu.CacheLinePad
	write     atomic.Int64
	_         cpu.CacheLinePad
	watermark atomic.Int64
	_         cpu.CacheLinePad
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

	n = copy(p.bufs[0][int(p.write.Load()):int(p.watermark.Load())], buf[0])

	wl := int64(len(buf[0]))
	if int64(len(p.bufs))-p.write.Load() >= wl {
		p.watermark.Store(p.write.Load() + wl)
		p.write.Add(wl)
	} else {
		if p.read.Load() < wl {
			// TODO Wait. Lock the lock-free :D
		}
		p.watermark.Store(p.write.Load())
		p.write.Store(wl)
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

	w, r := int(p.write.Load()), int(p.read.Load())
	m := int(p.watermark.Load())
	if w >= r {
		n = copy(buf[0], p.bufs[0][r:w])
		p.read.Store(p.write.Load())
	} else {
		n = copy(buf[0], p.bufs[0][r:m])
		p.read.Store(0)
	}

	return
}

var _ SignalReader = &altpipe{}
var _ SignalWriter = &altpipe{}

func Altpipe(nch int) (SignalReader, SignalWriter) {
	p := &altpipe{
		nch: nch,
	}

	return p, p
}
