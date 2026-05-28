package main

import (
	"os"
	"sync"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/waveform"
)

func experiments(n int, inputfile *os.File, output string) {
	if n == 1 {
		// Readall, bypass
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		wsw, err := wsr.MakeWriter(err, of)
		rr := dspio.TeeReader(wsr, wsw)
		d, err := dspio.ReadAll(err, rr)
		if err != nil {
			panic(err)
		}
		err = waveform.Dump(err, d[0])
		if err != nil {
			panic(err)
		}
	}
	if n == 2 {
		// Novelty curve
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		wsw, err := wsr.MakeWriter(err, of)
		if err != nil {
			panic(err)
		}
		dt := nanowarp.DetectorNew(1024, int(wsw.Format.SampleRate), 30, 300)
		err = dt.NoveltyCurveProcess(wsr, wsw)
		if err != nil {
			panic(err)
		}
	}
	if n == 3 {
		// Grain delay
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		wsw, err := wsr.MakeWriter(err, of)
		if err != nil {
			panic(err)
		}
		gr := dspio.NewGrainReader(1024, 1024, wsr)
		gw := dspio.NewGrainWriter(1024, 1024, wsw)
		buf := make([][]float64, wsr.NchRead())
		for ch := range buf {
			buf[ch] = make([]float64, 8192)
		}
		_, err = dspio.Copy(nil, gr, gw, buf)
		if err != nil {
			panic(err)
		}
	}
	if n == 4 {
		// The lack of grain delay in offline tract
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		wsw, err := wsr.MakeWriter(err, of)
		if err != nil {
			panic(err)
		}
		gr := dspio.NewOfflineGrainReader(1024, 1024, wsr)
		gw := dspio.NewOfflineGrainWriter(1024, 1024, wsw)
		buf := make([][]float64, wsr.NchRead())
		for ch := range buf {
			buf[ch] = make([]float64, 8192)
		}
		_, err = dspio.Copy(nil, gr, gw, buf)
		if err != nil {
			panic(err)
		}
	}
	if n == 5 {
		// Center-dilated novelty curve
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		wsw, err := wsr.MakeWriter(err, of)
		if err != nil {
			panic(err)
		}
		dt := nanowarp.DetectorNew(1024, int(wsw.Format.SampleRate), 30, 300)
		po, pi := dspio.GoPipe(2)
		wg := sync.WaitGroup{}
		wg.Add(2)
		go func() {
			defer wg.Done()
			// This piece of shit must be deferred after `defer wg.Done()` or else we get deadlock.
			// No, wg.Go() won't work.
			// Closing the pipe inside processors after EOF won't work either.
			defer pi.Close()
			_ = dt.NoveltyCurveProcess(wsr, pi)
		}()
		go func() {
			defer wg.Done()
			_ = dt.DilatePeakSelectProcess(po, wsw, 1, nil)
		}()
		wg.Wait()
	}
}
