package wav

import (
	"encoding/binary"
	"errors"
	"io"
	"unsafe"
)

var NotAWav = errors.New(`not a WAV file`)
var Malformed = errors.New(`malformed WAV file`)
var UnknownFormat = errors.New(`file is a WAV, but not in PCM or IEEE float format`)

type Format int

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
	WFormatTag          uint16
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

type Reader struct {
	dataseek   int64 // Start of sample data in bytes
	riff, data riffHeader
	fmtchunk

	// A list of all second-level RIFF headers identified in a file, except
	// `WAVE`, `fmt ` and `data`.
	Headermap []Cue
}

func (r *Reader) Properties() Properties {
	return Properties{
		Samples:    int(r.data.cksize) / int(r.NBlockAlign),
		Bytes:      int(r.riff.cksize) + 8, // Count RIFFxxxx 8-byte header too
		Samplerate: int(r.fmtchunk.NSamplesPerSec),
	}
}

// Cue is a cue point of a data found from its RIFF section in a file
// with byte index at which it was found.
//
// It points to the start of data.
type Cue struct {
	FourCC [4]byte
	Seek   int64
	Size   uint32
}

type riffHeader struct {
	fourcc [4]byte
	cksize uint32
}

func NewReader(rs io.ReadSeeker) (*Reader, error) {
	r := Reader{}

	var RIFF riffHeader
	err := binary.Read(rs, binary.LittleEndian, RIFF)
	if err != nil {
		return nil, err
	}
	if string(RIFF.fourcc[:]) != `RIFF` {
		return nil, NotAWav
	}

	var WAVE [4]byte
	err = binary.Read(rs, binary.LittleEndian, WAVE)
	if err != nil {
		return nil, err
	}
	if string(WAVE[:]) != `WAVE` {
		return nil, NotAWav
	}

	var fmt riffHeader
	err = binary.Read(rs, binary.LittleEndian, fmt)
	if err != nil {
		return nil, err
	}
	if string(fmt.fourcc[:]) != `fmt ` {
		return nil, Malformed
	}
	if !(fmt.cksize == 40 || fmt.cksize == 16 || fmt.cksize == 18) {
		return nil, Malformed
	}

	// Get seek before fmt chunk
	r.dataseek, err = rs.Seek(0, io.SeekCurrent)
	if err != nil {
		return nil, err
	}
	err = binary.Read(rs, binary.LittleEndian, r.fmtchunk)
	if err != nil {
		return nil, err
	}
	// Limit fmt chunk data to its real size.
	const fmtchunklen = unsafe.Sizeof(r.fmtchunk)
	clear((*[fmtchunklen]byte)(unsafe.Pointer(&r.fmtchunk))[:fmtchunklen-uintptr(fmt.cksize)])
	// Skip fmt chunk correctly.
	r.dataseek, err = rs.Seek(r.dataseek+int64(fmt.cksize), io.SeekStart)
	if err != nil {
		return nil, err
	}

	for {
		var data riffHeader
		err = binary.Read(rs, binary.LittleEndian, data)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if string(fmt.fourcc[:]) == `data` {

		}
		break
	}

	return &r, nil
}
