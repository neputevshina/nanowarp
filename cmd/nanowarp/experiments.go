package main

import (
	"os"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/waveform"
)

func experiments(n int, inputfile *os.File, output string) {
	if n == 1 {
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		d, err := dspio.ReadAll(err, wsr)
		if err != nil {
			panic(err)
		}
		err = waveform.Dump(err, d[0])
		if err != nil {
			panic(err)
		}
	}
	if n == 2 {
		of, err := os.Create(output)
		defer of.Close()
		wsr, err := NewWavSignalReader(err, inputfile)
		wsw, err := wsr.MakeWriter(err, of)
		if err != nil {
			panic(err)
		}
		dt := nanowarp.DetectorNew(512, 48000, 30, 300)
		err = dt.OnsetFunctionWriter(wsr, wsw)
		if err != nil {
			panic(err)
		}
	}
}
