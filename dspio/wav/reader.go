package wav

import (
	"encoding/binary"
	"io"

	"github.com/neputevshina/nanowarp/dspio"
	"github.com/zaf/g711"
	"golang.org/x/exp/constraints"
)

// Reader is a reader object for WAV files.
type Reader struct {
	readbuf []byte
	rs      io.ReadSeeker

	data       Cue // Start of sample data in bytes
	riff       riffHeader
	bytespersa int // Bytes per each single-channel sample
	fmtchunk   fmtchunk

	// A list of all second-level RIFF headers identified in a file,
	// except `WAVE` and `fmt `.
	Headermap []Cue
}

// Properties returns all relevant properties of a file.
func (r *Reader) Properties() Properties {
	return Properties{
		Samples:    int(r.data.Size) / int(r.fmtchunk.NBlockAlign),
		Bytes:      int(r.riff.Cksize) + 8, // Count RIFFxxxx 8-byte header too
		Samplerate: int(r.fmtchunk.NSamplesPerSec),
		Nch:        int(r.fmtchunk.NChannels),
		Format:     r.fmtchunk.WFormatTag,
	}
}

// NewReader creates a new Reader.
func NewReader(rs io.ReadSeeker) (*Reader, error) {
	r := Reader{}
	r.rs = rs

	err := binary.Read(rs, binary.LittleEndian, &r.riff)
	if err != nil {
		return nil, err
	}
	if string(r.riff.Fourcc[:]) != `RIFF` {
		return nil, NotAWav
	}

	var WAVE [4]byte
	err = binary.Read(rs, binary.LittleEndian, &WAVE)
	if err != nil {
		return nil, err
	}
	if string(WAVE[:]) != `WAVE` {
		return nil, NotAWav
	}

	var fmt riffHeader
	err = binary.Read(rs, binary.LittleEndian, &fmt)
	if err != nil {
		return nil, err
	}
	if string(fmt.Fourcc[:]) != `fmt ` {
		return nil, Malformed
	}

	// Get seek before fmt chunk.
	pre, err := rs.Seek(0, io.SeekCurrent)
	if err != nil {
		return nil, err
	}
	err = binary.Read(rs, binary.LittleEndian, &r.fmtchunk)
	if err != nil {
		return nil, err
	}
	// Cut unpopulated data.
	if fmt.Cksize == 16 || fmt.Cksize == 18 {
		r.fmtchunk.CbSize = 0
		r.fmtchunk.WValidBitsPerSample = 0
		r.fmtchunk.DwChannelMask = 0
		r.fmtchunk.SubFormat = [16]byte{}
	} else if fmt.Cksize != 40 {
		// Only allowed sizes for fmt chunk are 16, 18 and 40 bytes.
		return nil, Malformed
	}

	// Skip fmt chunk correctly.
	_, err = rs.Seek(pre+int64(fmt.Cksize), io.SeekStart)
	if err != nil {
		return nil, err
	}
	wavfmt := Format(r.fmtchunk.WFormatTag)
	if !(wavfmt == FormatPCM || wavfmt == FormatFloat) {
		return nil, UnsupportedFormat
	}
	if wavfmt == FormatFloat {
		if !(r.fmtchunk.WBitsPerSample == 32 || r.fmtchunk.WBitsPerSample == 64) {
			return nil, Malformed
		}
	}
	r.bytespersa = int(r.fmtchunk.WBitsPerSample >> 3)
	if r.fmtchunk.WBitsPerSample&7 > 0 {
		r.bytespersa++
	}

	for {
		var data riffHeader
		err = binary.Read(rs, binary.LittleEndian, &data)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		seek, err := rs.Seek(0, io.SeekCurrent)
		if err != nil {
			return nil, err
		}

		q := Cue{
			Seek:   seek,
			FourCC: data.Fourcc,
			Size:   data.Cksize,
		}
		r.Headermap = append(r.Headermap, q)
		if string(data.Fourcc[:]) == `data` {
			r.data = q
		}

		_, err = rs.Seek(int64(data.Cksize), io.SeekCurrent)
		if err != nil {
			if r.data.Seek != 0 && err == io.EOF {
				// Malformed, but workable.
				break
			} else {
				return nil, err
			}
		}
	}

	_, err = rs.Seek(r.data.Seek, io.SeekStart)
	if err != nil {
		return nil, err
	}
	return &r, nil
}

var _ dspio.SignalReader = &Reader{}

// NchRead returns the number of channels in a WAV file.
func (r *Reader) NchRead() int {
	return int(r.fmtchunk.NChannels)
}

// SignalRead reads a non-overlapped multichannel grain of size len(buf[0]) from the input.
func (r *Reader) SignalRead(prr error, buf [][]float64) (n int, err error) {
	blen := len(buf[0]) * int(r.fmtchunk.NBlockAlign)
	if len(r.readbuf) < blen {
		r.readbuf = make([]byte, blen)
	}
	r.readbuf = r.readbuf[:blen]

	nbs, err := r.rs.Read(r.readbuf)
	n = nbs / int(r.fmtchunk.NBlockAlign)
	if err != nil {
		return n, err
	}

	switch r.fmtchunk.WFormatTag {
	case FormatALaw:
		decodeCompanded(r, r.readbuf, buf, 0)
	case FormatMuLaw:
		decodeCompanded(r, r.readbuf, buf, 1)
	case FormatPCM:
		decodeInteger(r, int(r.fmtchunk.WBitsPerSample), r.readbuf, buf)
	case FormatFloat:
		switch r.fmtchunk.WBitsPerSample {
		case 32:
			decodeFloat[float32](r, r.readbuf, buf)
		case 64:
			decodeFloat[float64](r, r.readbuf, buf)
		default:
			panic(`fatal: early WAV file format verification failed`)
		}
	}
	return n, nil
}

func decodeCompanded(r *Reader, bbuf []byte, buf [][]float64, law int) {
	var sa byte
	nch := int(r.fmtchunk.NChannels)
	for i := range buf[0] {
		for ch := range nch {
			sl := bbuf[(i*nch+ch)*int(r.fmtchunk.NBlockAlign)/r.bytespersa:]
			binary.Decode(sl, binary.LittleEndian, &sa)
			dec := []func(byte) int16{g711.DecodeAlawFrame, g711.DecodeUlawFrame}[law]
			buf[ch][i] = float64(dec(sa)) / float64(1<<15)
		}
	}
}

func decodeInteger(r *Reader, nbits int, bbuf []byte, buf [][]float64) {
	var sasa [8]byte
	var sa uint64
	nch := int(r.fmtchunk.NChannels)
	for i := range buf[0] {
		for ch := range nch {
			copy(sasa[:], bbuf[(i*nch+ch)*int(r.fmtchunk.NBlockAlign)/r.bytespersa:])
			binary.Decode(sasa[:], binary.LittleEndian, &sa)
			sa += 1 << (nbits - 1)
			sa &= uint64(1<<nbits - 1)
			buf[ch][i] = float64(sa)/float64(int(1)<<(nbits-1)) - 1
		}
	}
}

func decodeFloat[T constraints.Float](r *Reader, bbuf []byte, buf [][]float64) {
	var sa T
	for i := range buf {
		for ch := range int(r.fmtchunk.NChannels) {
			binary.Decode(bbuf[i+ch*int(r.fmtchunk.NChannels):][:int(r.fmtchunk.NBlockAlign)/int(r.fmtchunk.NChannels)],
				binary.LittleEndian, &sa)
			buf[ch][i] = float64(sa)
		}
	}
}
