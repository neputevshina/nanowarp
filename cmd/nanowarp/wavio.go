package main

import (
	"github.com/neputevshina/nanowarp/dspio"
	"github.com/neputevshina/nanowarp/wav"
	"github.com/youpy/go-riff"
)

type WavSignalReader struct {
	*wav.Reader
	*wav.WavFormat

	buf []wav.Sample
}

var _ dspio.SignalReader = &WavSignalReader{}

func NewWavSignalReader(prr error, file riff.RIFFReader) (wsr *WavSignalReader, err error) {
	if prr != nil {
		return nil, prr
	}
	wsr = &WavSignalReader{}
	wsr.Reader = wav.NewReader(file)
	wsr.WavFormat, err = wsr.Reader.Format()
	if err != nil {
		return nil, err
	}
	return
}

func (w *WavSignalReader) NchRead() int {
	return int(w.WavFormat.NumChannels)
}

func (w *WavSignalReader) SignalRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}

	if len(w.buf) < len(buf[0]) {
		w.buf = make([]wav.Sample, len(buf[0]))
	}
	w.buf = w.buf[:len(buf[0])]
	n, _, err = w.Reader.ReadSamples(w.buf)
	for i := range n {
		for ch := range w.NumChannels {
			buf[ch][i] = float64(w.buf[i].Values[ch])
		}
	}
	return
}

type WavSignalWriter struct{}

func (w *WavSignalWriter) NchWrite() int {
	return 0
}

func (w *WavSignalWriter) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	return
}

var _ dspio.SignalWriter = &WavSignalWriter{}
