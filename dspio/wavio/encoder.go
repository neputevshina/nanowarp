package wavio

import (
	"encoding/binary"
	"io"
	"slices"
	"unsafe"

	"github.com/neputevshina/nanowarp/dspio"
)

// Encoder is a writer of a WAV file.
//
// The only supported format is IEEE float32.
type Encoder struct {
	wr      io.WriteSeeker
	fs, nch int

	blockalign   int
	dataseek     int64 // `data` section seek
	riffc, datac int64 // Samples written so far

	bytebuf []byte
}

// NewEncoder creates a new IEEE 32-bit float WAV file encoder with specified format
// and writes the WAV header to file.
//
// The file is assumed to be empty or reset to seek 0 before writing.
//
// format must be [FormatFloat] and bits must be 32.
// In future, other formats may be supported.
func NewEncoder(file io.WriteSeeker, samplerate int, channels int, format Format, bits int) (e *Encoder, err error) {
	e = &Encoder{
		wr:  file,
		fs:  samplerate,
		nch: channels,
	}

	if !(format == FormatFloat && bits == 32) {
		panic(`incorrect format: only 32-bit float supported`)
	}

	// TODO
	switch format {
	case FormatPCM:
		if bits < 1 && bits > 32 {
			panic(`incorrect PCM format: number of bits must be between 1 and 32`)
		}
	case FormatFloat:
		if !(bits == 32 || bits == 64) {
			panic(`incorrect float format: number of bits must be either 32 or 64`)
		}
	default:
		panic(`can't write this format`)
	}

	e.blockalign = (bits>>3 + boolint(bits&7 > 0)) * channels
	bytespersec := samplerate * e.blockalign
	if int(uint32(bytespersec)) != bytespersec {
		panic(`incorrect format: >4 GiB per second, which is impossible`)
	}

	// Last four bytes are uint32(10) in little-endian, the size of a format chunk.
	prelude := []byte("RIFF\x00\x00\x00\x00WAVEfmt \x10\x00\x00\x00")
	n, err := file.Write(prelude)
	if err != nil {
		return nil, err
	} else if n != len(prelude) {
		return nil, io.ErrShortWrite
	}

	fmt := fmtshortchunk{
		WFormatTag:      format,
		NChannels:       uint16(channels),
		NSamplesPerSec:  uint32(samplerate),
		NAvgBytesPerSec: uint32(bytespersec),
		NBlockAlign:     uint16(e.blockalign),
		WBitsPerSample:  uint16(bits),
	}
	err = binary.Write(file, binary.LittleEndian, fmt)
	if err != nil {
		return nil, err
	}

	e.riffc = int64(n) + int64(unsafe.Sizeof(fmt))
	return
}

// WriteRiffChunk copies RIFF chunk data from src until it encounters an EOF.
// Resulting section size will be set to total number of bytes read from src.
//
// WriteRiffChunk can't be called after the first call to [Encoder.SignalWrite].
//
// If total number of bytes written is greater than 4 GiB,
// WriteRiffChunk returns [ErrOverflow].
func (e *Encoder) WriteRiffChunk(fourcc [4]byte, src io.Reader) error {
	if e.dataseek != 0 {
		panic(`can't write RIFF chunk after sample data has started being written`)
	}

	n, err := e.wr.Write(fourcc[:])
	if err != nil {
		return err
	} else if n != 4 {
		return io.ErrShortWrite
	}
	sizeseek, err := e.wr.Seek(0, io.SeekCurrent)
	if err != nil {
		return err
	}

	n, err = e.wr.Write((&[4]byte{})[:])
	if err != nil {
		return err
	} else if n != 4 {
		return io.ErrShortWrite
	}

	nn, err := io.Copy(e.wr, io.LimitReader(src, 1<<32))
	e.riffc += nn + 8
	if err != nil {
		return err
	}
	if e.riffc > 1<<32 {
		return ErrOverflow
	}

	nextseek, err := e.wr.Seek(0, io.SeekCurrent)
	if err != nil {
		return err
	}

	_, err = e.wr.Seek(sizeseek, io.SeekStart)
	if err != nil {
		return err
	}
	err = binary.Write(e.wr, binary.LittleEndian, uint32(nn))
	if err != nil {
		return err
	}

	// All byte seek addresses in WAV must be even.
	// If the size of a section is odd, we still write the real size of a
	// section, but advance the seek to even address in bytes.
	_, err = e.wr.Seek(even(nextseek), io.SeekStart)
	return err
}

// TODO
//
// WriteInfo writes LIST INFO section, storing the information about music
// inside the file.
// Despite it accepting a map from FourCC to multiple sections, it writes only
// one section per each FourCC, concatenating each chunk of data with \n, because
// having more than one subsection of each type in INFO is non-standard.
//
// WriteInfo can't be called after the first call to [Encoder.SignalWrite].
//
// If total number of bytes written is greater than 4 GiB,
// WriteInfo returns [ErrOverflow].
func (e *Encoder) WriteInfo(data map[[4]byte][]string) error {
	panic(`unimplemented`)
}

var _ dspio.SignalWriter = &Encoder{}

// NchWrite returns the number of channels this encoder will write.
func (e *Encoder) NchWrite() int { return e.nch }

// SignalWrite writes a non-overlapping multichannel grain from buf.
//
// If total number of bytes written is greater than 4 GiB, SignalWrite returns
// [ErrOverflow], but writes as many samples as it can and forcefully closes
// the Encoder.
//
// In the case of SignalWrite, [ErrOverflow] must be treated similar to EOF.
func (e *Encoder) SignalWrite(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}

	if e.dataseek == 0 {
		e.dataseek, err = e.wr.Seek(0, io.SeekCurrent)
		if err != nil {
			return 0, err
		}
		n, err := e.wr.Write([]byte("data\x00\x00\x00\x00"))
		if err != nil {
			return 0, err
		} else if n != 8 {
			return 0, io.ErrShortWrite
		}

	}

	nb := len(buf[0]) * e.blockalign
	if e.riffc+int64(nb) > 1<<32 {
		return 0, ErrOverflow
	}
	if len(e.bytebuf) != nb {
		e.bytebuf = slices.Grow(e.bytebuf[:0], len(buf)*len(buf[0])*e.blockalign)[:nb]
	}

	pn := 0
	for j := range buf[0] {
		for ch := range buf {
			n, err = binary.Encode(e.bytebuf[pn:], binary.LittleEndian, float32(buf[ch][j]))
			if err != nil {
				panic(`unreachable`)
			}
			pn += n
		}
	}

	n, err = e.wr.Write(e.bytebuf)
	if n != len(e.bytebuf) {
		err = io.ErrShortWrite
	}

	e.riffc += int64(n)
	e.datac += int64(n)

	return
}

// Close finalizes the WAV file write, setting the total size of file and
// sample data sections in headers to amount of data gone through Encoder.
//
// It does not attempt to close the underlying [io.WriteSeeker].
func (e *Encoder) Close() (err error) {
	_, err = e.wr.Seek(4, io.SeekStart)
	if err != nil {
		return
	}
	err = binary.Write(e.wr, binary.LittleEndian, uint32(e.riffc))
	if err != nil {
		return
	}

	_, err = e.wr.Seek(e.dataseek+4, io.SeekStart)
	if err != nil {
		return
	}
	err = binary.Write(e.wr, binary.LittleEndian, uint32(e.datac))
	if err != nil {
		return
	}

	_, err = e.wr.Seek(e.datac, io.SeekCurrent)
	return
}

func boolint(x bool) int {
	if x {
		return 1
	}
	return 0
}
