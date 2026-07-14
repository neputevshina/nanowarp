package dspio

// TODO Write own WAV codec, without all of the peculiarities of youpy's code.

import (
	"fmt"
	"io"
	"math"

	"github.com/neputevshina/nanowarp/wav"
	"github.com/youpy/go-riff"
)

type WavSignalReader struct {
	*wav.Reader
	*wav.WavFormat

	buf []wav.Sample
}

var _ SignalReader = &WavSignalReader{}

func NewWavSignalReader(prr error, file riff.RIFFReader) (wsr *WavSignalReader, err error) {
	wsr = &WavSignalReader{}
	if prr != nil {
		return nil, prr
	}
	wsr.Reader = wav.NewReader(file)
	_, _ = wsr.Duration() // WavData is populated by an unexported function.
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
			switch w.WavFormat.AudioFormat {
			case wav.AudioFormatIEEEFloat:
				switch w.BitsPerSample {
				case 32:
					buf[ch][i] = float64(math.Float32frombits(uint32(w.buf[i].Values[ch])))
				case 64:
					buf[ch][i] = float64(math.Float64frombits(uint64(w.buf[i].Values[ch])))
				default:
					return n, fmt.Errorf("malformed wav: float bit depth is not 32 nor 64 (got %v)", w.BitsPerSample)
				}
			default:
				norm := math.Pow(2, float64(w.BitsPerSample-1))
				buf[ch][i] = float64(w.buf[i].Values[ch]) / norm
			}

		}
	}
	return
}

func (w *WavSignalReader) MakeWriter(prr error, of io.Writer) (wsw *WavSignalWriter, err error) {
	if prr != nil {
		return nil, prr
	}
	return NewWavSignalWriter(err, of, int(w.Size)/int(w.NumChannels)/int(w.BitsPerSample/8), int(w.NumChannels), int(w.SampleRate))
}

type WavSignalWriter struct {
	*wav.Writer
}

var _ SignalWriter = &WavSignalWriter{}

func NewWavSignalWriter(prr error, file io.Writer, length int, nch int, fs int) (wsw *WavSignalWriter, err error) {
	if prr != nil {
		return nil, prr
	}
	wsw = &WavSignalWriter{}
	wsw.Writer = wav.NewWriter(file, uint32(length), uint16(nch), uint32(fs), 32, true)
	return
}

func (w *WavSignalWriter) NchWrite() int {
	return int(w.Writer.Format.NumChannels)
}

func (w *WavSignalWriter) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}

	nbuf := len(buf[0])
	sabuf := make([]wav.Sample, 0, nbuf)

	for i := range nbuf {
		sa := wav.Sample{}
		for ch := range w.NchWrite() {
			sa.Values[ch] = int(math.Float32bits(float32(buf[ch][i])))
		}
		sabuf = append(sabuf, sa)
	}
	err = w.WriteSamples(sabuf)
	if err != nil {
		panic(err)
	}
	n = len(sabuf)

	return
}
