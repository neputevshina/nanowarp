package wavio

import (
	"io"

	"github.com/neputevshina/nanowarp/dspio"
)

// Encoder is a writer of a WAV file.
//
// The only supported format is IEEE float32.
type Encoder struct {
	wr      io.WriteSeeker
	fs, nch int
}

// NewEncoder creates a new IEEE 32-bit float WAV file encoder with specified sample rate and
// number of channels and writes the WAV header to file.
func NewEncoder(file io.WriteSeeker, samplerate int, channels int) (*Encoder, error) {
	panic(`unimplemented`)
	// e := Encoder{
	// 	wr:  file,
	// 	fs:  samplerate,
	// 	nch: channels,
	// }
	// return &e, nil
}

// Close finalizes the WAV file write, setting the total size of file and
// sample data sections to amount of data gone through Encoder.
//
// Underlying [io.WriteSeeker] is not closed.
func (e *Encoder) Close() error {
	panic(`unimplemented`)
}

// WriteRiffChunk writes a RIFF chunk data from src.
//
// It copies data from the src until it encounters an EOF.
// Resulting section size will be set to total number of bytes read from src.
//
// WriteRiffChunk can't be called after signal writing has started.
//
// If total number of bytes written is greater than 4 GiB,
// WriteRiffChunk returns [ErrOverflow].
func (e *Encoder) WriteRiffChunk(fourcc [4]byte, src io.Reader) error {
	panic(`unimplemented`)
}

var _ dspio.SignalWriter = &Encoder{}

// NchWrite returns the number of channels this encoder will write.
func (e *Encoder) NchWrite() int {
	panic("unimplemented")
}

// SignalWrite writes a non-overlapping multichannel grain from buf.
//
// If total number of bytes written is greater than 4 GiB, SignalWrite returns
// [ErrOverflow], but writes as many samples as it can and forcefully closes
// the Encoder.
//
// In the case of SignalWrite, [ErrOverflow] must be treated similar to EOF.
func (e *Encoder) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	panic("unimplemented")
}
