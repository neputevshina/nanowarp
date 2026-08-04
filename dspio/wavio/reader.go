package wavio

import (
	"encoding/binary"
	"io"

	"github.com/neputevshina/nanowarp/dspio"
	"github.com/zaf/g711"
	"golang.org/x/exp/constraints"
)

// Decoder is a reader object for WAV files.
type Decoder struct {
	readbuf []byte
	rs      io.ReadSeeker

	data       Section // Start of sample data in bytes
	riff       riffHeader
	bytespersa int // Bytes per each single-channel sample
	fmtchunk   fmtchunk

	// List of all second-level RIFF headers identified in a file.
	Headermap []Section
}

// Properties returns all relevant properties of a file.
func (r *Decoder) Properties() Properties {
	return Properties{
		Samples:    int(r.data.Size) / int(r.fmtchunk.NBlockAlign),
		Bytes:      int(r.riff.Cksize) + 8, // Count RIFFxxxx 8-byte header too
		Samplerate: int(r.fmtchunk.NSamplesPerSec),
		Nch:        int(r.fmtchunk.NChannels),
		Format:     r.fmtchunk.WFormatTag,
	}
}

// NewDecoder creates a new Reader.
func NewDecoder(rs io.ReadSeeker) (*Decoder, error) {
	r := Decoder{}
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
	q := Section{
		Seek:   pre,
		FourCC: fmt.Fourcc,
		Size:   fmt.Cksize,
	}
	r.Headermap = append(r.Headermap, q)

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
	// Validate format.
	wavfmt := Format(r.fmtchunk.WFormatTag)
	if !(wavfmt == FormatPCM || wavfmt == FormatFloat ||
		wavfmt == FormatALaw || wavfmt == FormatMuLaw) {
		return nil, UnsupportedFormat
	}
	if wavfmt == FormatFloat {
		if !(r.fmtchunk.WBitsPerSample == 32 || r.fmtchunk.WBitsPerSample == 64) {
			return nil, Malformed
		}
	}
	if wavfmt == FormatALaw || wavfmt == FormatMuLaw {
		if r.fmtchunk.WBitsPerSample != 8 {
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

		q := Section{
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

// Rewind seeks the file to the start of sample data.
//
// If you need the seek index itself, search for `data` section in r.Headermap.
func (r *Decoder) Rewind() (err error) {
	_, err = r.rs.Seek(r.data.Seek, io.SeekStart)
	return
}

var _ dspio.SignalReader = &Decoder{}

// NchRead returns the number of channels in a WAV file.
func (r *Decoder) NchRead() int {
	return int(r.fmtchunk.NChannels)
}

// SignalRead reads a non-overlapped multichannel grain of size len(buf[0]) from the input.
func (r *Decoder) SignalRead(prr error, buf [][]float64) (n int, err error) {
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
		decodeCompanded(r, r.readbuf, buf, false)
	case FormatMuLaw:
		decodeCompanded(r, r.readbuf, buf, true)
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

func decodeCompanded(r *Decoder, bbuf []byte, buf [][]float64, ulaw bool) {
	var sa byte
	nch := int(r.fmtchunk.NChannels)
	for i := range buf[0] {
		for ch := range nch {
			sa = bbuf[(i*nch+ch)*int(r.fmtchunk.NBlockAlign)/r.bytespersa:][0]
			if ulaw {
				buf[ch][i] = float64(g711.DecodeUlawFrame(sa)) / float64(1<<15)
			} else {
				buf[ch][i] = float64(g711.DecodeAlawFrame(sa)) / float64(1<<15)
			}

		}
	}
}

func decodeInteger(r *Decoder, nbits int, bbuf []byte, buf [][]float64) {
	var sasa [8]byte
	var sa uint64
	nch := int(r.fmtchunk.NChannels)
	for i := range buf[0] {
		for ch := range nch {
			copy(sasa[:], bbuf[(i*nch+ch)*int(r.fmtchunk.NBlockAlign)/nch:])
			// NOTE(neputevshina): This call is inlined.
			binary.Decode(sasa[:], binary.LittleEndian, &sa)
			sa += 1 << (nbits - 1)
			sa &= uint64(1<<nbits - 1)
			buf[ch][i] = float64(sa)/float64(int(1)<<(nbits-1)) - 1
		}
	}
}

func decodeFloat[T constraints.Float](r *Decoder, bbuf []byte, buf [][]float64) {
	var sa T
	nch := int(r.fmtchunk.NChannels)
	for i := range buf[0] {
		for ch := range int(r.fmtchunk.NChannels) {
			// NOTE(neputevshina): This call is inlined.
			binary.Decode(bbuf[(i*nch+ch)*int(r.fmtchunk.NBlockAlign)/nch:],
				binary.LittleEndian, &sa)
			buf[ch][i] = float64(sa)
		}
	}
}

// InfoChunk reads the INFO chunk from the file if there is one.
// If chunk is not found, data is nil.
//
// InfoChunk assumes that there is none or only one INFO chunk in a file, but could be many chunks inside of it.
//
// If finished without error, it automatically rewinds r after reading the chunk.
//
// More advanced user might consider reading ID3 tags from files.
// Those can be read manually using chunk pointer information in r.Headermap.
func (r *Decoder) InfoChunk() (data map[[4]byte][]string, err error) {
	siz := int64(0)
	for _, m := range r.Headermap {
		siz = int64(m.Size)
		if string(m.FourCC[:]) == `LIST` {
			_, err := r.rs.Seek(m.Seek, io.SeekStart)
			if err != nil {
				return nil, err
			}

			var sub [4]byte
			err = binary.Read(r.rs, binary.LittleEndian, &sub)
			if err != nil {
				return nil, err
			}
			if string(sub[:]) == `INFO` {
				goto found
			}
		}
	}
	return nil, r.Rewind()

found:
	data = make(map[[4]byte][]string)
	for toll := int64(siz) - 4; toll > 0; {
		var head riffHeader
		err := binary.Read(r.rs, binary.LittleEndian, &head)
		if err != nil {
			return nil, err
		}
		toll -= 8
		str := make([]byte, even(head.Cksize))
		_, err = r.rs.Read(str)
		if err != nil {
			return nil, err
		}
		data[head.Fourcc] = append(data[head.Fourcc], string(str[:head.Cksize]))
		toll -= int64(even(head.Cksize))
	}

	return data, r.Rewind()
}

func even[T constraints.Integer](x T) T {
	return x + x%2
}
