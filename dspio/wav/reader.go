package wav

import (
	"encoding/binary"
	"io"
)

type Reader struct {
	data Cue // Start of sample data in bytes
	riff riffHeader
	fmtchunk

	// A list of all second-level RIFF headers identified in a file,
	// except `WAVE` and `fmt `.
	Headermap []Cue
}

func (r *Reader) Properties() Properties {
	return Properties{
		Samples:    int(r.data.Size) / int(r.NBlockAlign),
		Bytes:      int(r.riff.cksize) + 8, // Count RIFFxxxx 8-byte header too
		Samplerate: int(r.fmtchunk.NSamplesPerSec),
	}
}

func NewReader(rs io.ReadSeeker) (*Reader, error) {
	r := Reader{}

	err := binary.Read(rs, binary.LittleEndian, r.riff)
	if err != nil {
		return nil, err
	}
	if string(r.riff.fourcc[:]) != `RIFF` {
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

	// Get seek before fmt chunk.
	pre, err := rs.Seek(0, io.SeekCurrent)
	if err != nil {
		return nil, err
	}
	err = binary.Read(rs, binary.LittleEndian, r.fmtchunk)
	if err != nil {
		return nil, err
	}
	// Cut unpopulated data.
	if fmt.cksize == 16 || fmt.cksize == 18 {
		r.fmtchunk.CbSize = 0
		r.fmtchunk.WValidBitsPerSample = 0
		r.fmtchunk.DwChannelMask = 0
		r.fmtchunk.SubFormat = [16]byte{}
	} else if fmt.cksize != 40 {
		// Only allowed sizes for fmt chunk are 16, 18 and 40 bytes.
		return nil, Malformed
	}

	// Skip fmt chunk correctly.
	_, err = rs.Seek(pre+int64(fmt.cksize), io.SeekStart)
	if err != nil {
		return nil, err
	}
	wavfmt := Format(r.fmtchunk.WFormatTag)
	if !(wavfmt == FormatPCM || wavfmt == FormatFloat) {
		return nil, UnsupportedFormat
	}
	if wavfmt == FormatFloat {
		if !(r.WBitsPerSample == 32 || r.WBitsPerSample == 64) {
			return nil, Malformed
		}
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
		seek, err := rs.Seek(0, io.SeekCurrent)
		if err != nil {
			return nil, err
		}

		q := Cue{
			Seek:   seek,
			FourCC: data.fourcc,
			Size:   data.cksize,
		}
		r.Headermap = append(r.Headermap, q)
		if string(data.fourcc[:]) == `data` {
			r.data = q
		}

		_, err = rs.Seek(int64(data.cksize), io.SeekCurrent)
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
