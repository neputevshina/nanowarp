// Interpolated from Go standard library source.
// src/io/pipe.go

// TODO Replace with lock-free SPSC queue.

// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package dspio

import (
	"fmt"
	"io"
	"sync"
)

// onceError is an object that will only store an error once.
type onceError struct {
	sync.Mutex // guards following
	err        error
}

func (a *onceError) Store(err error) {
	a.Lock()
	defer a.Unlock()
	if a.err != nil {
		return
	}
	a.err = err
}
func (a *onceError) Load() error {
	a.Lock()
	defer a.Unlock()
	return a.err
}

// A gopipe is the shared gopipe structure underlying PipeReader and PipeWriter.
type gopipe struct {
	wrMu sync.Mutex // Serializes Write operations
	wrCh chan [][]float64
	rdCh chan int
	nch  int

	once sync.Once // Protects closing done
	done chan struct{}
	rerr onceError
	werr onceError
}

func (p *gopipe) read(b [][]float64) (n int, err error) {
	select {
	case <-p.done:
		return 0, p.readCloseError()
	default:
	}

	select {
	case bw := <-p.wrCh:
		nr := 0
		for ch := 0; ch < p.nch; ch++ {
			nr = copy(b[ch], bw[ch])
		}
		p.rdCh <- nr
		return nr, nil
	case <-p.done:
		return 0, p.readCloseError()
	}
}

func (p *gopipe) closeRead(err error) error {
	if err == nil {
		err = io.ErrClosedPipe
	}
	p.rerr.Store(err)
	p.once.Do(func() { close(p.done) })
	return nil
}

func (p *gopipe) write(b [][]float64) (n int, err error) {
	select {
	case <-p.done:
		return 0, p.writeCloseError()
	default:
		p.wrMu.Lock()
		defer p.wrMu.Unlock()
	}

	for once := true; once || len(b) > 0; once = false {
		select {
		case p.wrCh <- b:
			nw := <-p.rdCh
			for ch := 0; ch < p.nch; ch++ {
				b[ch] = b[ch][nw:]
			}
			n += nw
		case <-p.done:
			return n, p.writeCloseError()
		}
	}
	return n, nil
}

func (p *gopipe) closeWrite(err error) error {
	if err == nil {
		err = io.EOF
	}
	p.werr.Store(err)
	p.once.Do(func() { close(p.done) })
	return nil
}

// readCloseError is considered internal to the pipe type.
func (p *gopipe) readCloseError() error {
	rerr := p.rerr.Load()
	if werr := p.werr.Load(); rerr == nil && werr != nil {
		return werr
	}
	return io.ErrClosedPipe
}

// writeCloseError is considered internal to the pipe type.
func (p *gopipe) writeCloseError() error {
	werr := p.werr.Load()
	if rerr := p.rerr.Load(); werr == nil && rerr != nil {
		return rerr
	}
	return io.ErrClosedPipe
}

// A goPipeReader is the read half of a pipe.
type goPipeReader struct{ gopipe }

// Read implements the standard Read interface:
// it reads data from the pipe, blocking until a writer
// arrives or the write end is closed.
// If the write end is closed with an error, that error is
// returned as err; otherwise err is EOF.
func (r *goPipeReader) Read(prr error, data [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(data) != r.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, r.nch, len(data)))
	}
	return r.gopipe.read(data)
}

// Close closes the reader; subsequent writes to the
// write half of the pipe will return the error [ErrClosedPipe].
func (r *goPipeReader) Close() error {
	return r.CloseWithError(nil)
}

// CloseWithError closes the reader; subsequent writes
// to the write half of the pipe will return the error err.
//
// CloseWithError never overwrites the previous error if it exists
// and always returns nil.
func (r *goPipeReader) CloseWithError(err error) error {
	return r.gopipe.closeRead(err)
}

// A goPipeWriter is the write half of a pipe.
type goPipeWriter struct{ r goPipeReader }

// Write implements the standard Write interface:
// it writes data to the pipe, blocking until one or more readers
// have consumed all the data or the read end is closed.
// If the read end is closed with an error, that err is
// returned as err; otherwise err is [ErrClosedPipe].
func (w *goPipeWriter) Write(prr error, data [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	if len(data) != w.r.nch {
		panic(fmt.Errorf(`different number of channels: expected %d, got %d`, w.r.nch, len(data)))
	}
	return w.r.gopipe.write(data)
}

// Close closes the writer; subsequent reads from the
// read half of the pipe will return no bytes and EOF.
func (w *goPipeWriter) Close() error {
	return w.CloseWithError(nil)
}

// CloseWithError closes the writer; subsequent reads from the
// read half of the pipe will return no bytes and the error err,
// or EOF if err is nil.
//
// CloseWithError never overwrites the previous error if it exists
// and always returns nil.
func (w *goPipeWriter) CloseWithError(err error) error {
	return w.r.gopipe.closeWrite(err)
}

func goPipe(channels int) (*goPipeReader, *goPipeWriter) {
	pw := &goPipeWriter{r: goPipeReader{gopipe: gopipe{
		wrCh: make(chan [][]float64),
		rdCh: make(chan int),
		nch:  channels,
		done: make(chan struct{}),
	}}}
	return &pw.r, pw
}
