package wav

import "errors"

var NotAWav = errors.New(`not a WAV file`)
var Malformed = errors.New(`malformed WAV file`)
var UnsupportedFormat = errors.New(`file is a WAV, but not in linear PCM or IEEE float format`)

type Format uint16

const (
	FormatPCM        Format = 0x0001
	FormatFloat             = 0x0003
	FormatALaw              = 0x0006
	FormatMuLaw             = 0x0007
	FormatExtensible        = 0xFFFE
)

type Properties struct {
	Samples    int    // Amount of multichannel samples in file
	Bytes      int    // Length of file in bytes per its RIFF header
	Samplerate int    // Samples per second
	Nch        int    // Number of channels
	Format     Format // Format of samples
}

type fmtchunk struct {
	WFormatTag          Format
	NChannels           uint16
	NSamplesPerSec      uint32
	NAvgBytesPerSec     uint32
	NBlockAlign         uint16
	WBitsPerSample      uint16
	CbSize              uint16
	WValidBitsPerSample uint16
	DwChannelMask       uint32
	SubFormat           [16]byte
}

// Cue is a pointer to a RIFF section in a file.
//
// Seek points to the start of data, start of the section itself is at Seek - 8.
type Cue struct {
	Seek   int64
	FourCC [4]byte
	Size   uint32
}

type riffHeader struct {
	fourcc [4]byte
	cksize uint32
}
