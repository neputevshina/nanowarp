package wavio

import (
	"errors"
)

var ErrNotAWav = errors.New(`not a WAV file`)
var ErrMalformed = errors.New(`malformed WAV file`)
var ErrUnsupportedFormat = errors.New(`file is a WAV, but not in linear PCM or IEEE float format`)

var ErrOverflow = errors.New(`total number of bytes written is greater than 4 GiB`)

type Format uint16

// Recognized and supported formats:
const (
	FormatPCM   Format = 0x0001
	FormatFloat        = 0x0003
	FormatALaw         = 0x0006
	FormatMuLaw        = 0x0007
)

// Properties contains all relevant WAV file properties.
type Properties struct {
	Samples    int    // Length of file in number of multichannel samples
	Bytes      int    // Length of file in bytes per its RIFF header
	Samplerate int    // Number of samples per second
	Nch        int    // Number of channels
	Format     Format // Format of numbers in file
}

type fmtchunk struct {
	fmtshortchunk
	CbSize              uint16
	WValidBitsPerSample uint16
	DwChannelMask       uint32
	SubFormat           [16]byte
}

type fmtshortchunk struct {
	WFormatTag      Format
	NChannels       uint16
	NSamplesPerSec  uint32
	NAvgBytesPerSec uint32
	NBlockAlign     uint16
	WBitsPerSample  uint16
}

// Section is a pointer to a RIFF section in a file.
//
// Seek points to the start of data, start of the section itself is 8 bytes before Seek.
type Section struct {
	Seek   int64
	FourCC [4]byte
	Size   uint32
}

type riffHeader struct {
	Fourcc [4]byte
	Cksize uint32
}
