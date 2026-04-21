package dspio

import (
	"errors"
	"io"
	"testing"
)

// mock implementations

type mockReader struct {
	data [][]float64
	pos  int
	hop  int
}

func (m *mockReader) NchRead() int { return len(m.data) }

func (m *mockReader) SignalRead(_ error, buf [][]float64) (int, error) {
	if m.pos+m.hop > len(m.data[0]) {
		return 0, io.EOF
	}
	for ch := range buf {
		copy(buf[ch], m.data[ch][m.pos:m.pos+m.hop])
	}
	m.pos += m.hop
	return m.hop, nil
}

type mockWriter struct {
	data [][]float64
}

func newMockWriter(nch int) *mockWriter {
	return &mockWriter{data: make([][]float64, nch)}
}

func (m *mockWriter) NchWrite() int { return len(m.data) }

func (m *mockWriter) SignalWrite(_ error, buf [][]float64) (int, error) {
	for ch := range buf {
		m.data[ch] = append(m.data[ch], buf[ch]...)
	}
	return len(buf[0]), nil
}

type errReader struct{ nch int }

func (e *errReader) NchRead() int { return e.nch }
func (e *errReader) SignalRead(_ error, _ [][]float64) (int, error) {
	return 0, errors.New("read error")
}

type errWriter struct{ nch int }

func (e *errWriter) NchWrite() int { return e.nch }
func (e *errWriter) SignalWrite(_ error, _ [][]float64) (int, error) {
	return 0, errors.New("write error")
}

// helpers

func makeGrainBuf(nch, n int) [][]float64 {
	g := make([][]float64, nch)
	for ch := range g {
		g[ch] = make([]float64, n)
	}
	return g
}

func slicesEqual(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// GrainReader

func TestGrainReaderPrr(t *testing.T) {
	r := NewGrainReader(4, 2, &mockReader{data: [][]float64{{1, 2, 3, 4}}, hop: 2})
	sentinel := errors.New("sentinel")
	n, err := r.SignalRead(sentinel, makeGrainBuf(1, 4))
	if n != 0 || err != sentinel {
		t.Fatalf("expected (0, sentinel), got (%d, %v)", n, err)
	}
}

func TestGrainReaderChannelPanic(t *testing.T) {
	r := NewGrainReader(4, 2, &mockReader{data: [][]float64{{1, 2, 3, 4}}, hop: 2})
	defer func() {
		if recover() == nil {
			t.Error("expected panic on channel mismatch")
		}
	}()
	r.SignalRead(nil, makeGrainBuf(2, 4)) // reader has nch=1, grain has 2
}

func TestGrainReaderNch(t *testing.T) {
	r := NewGrainReader(4, 2, &mockReader{data: [][]float64{{}, {}}, hop: 2})
	if r.NchRead() != 2 {
		t.Fatalf("expected 2, got %d", r.NchRead())
	}
}

func TestGrainReaderOverlap(t *testing.T) {
	// nfft=8, hop=4: each consecutive grain shares nfft-hop=4 samples
	data := make([]float64, 16)
	for i := range data {
		data[i] = float64(i + 1)
	}
	r := NewGrainReader(8, 4, &mockReader{data: [][]float64{data}, hop: 4})
	g1 := makeGrainBuf(1, 8)
	g2 := makeGrainBuf(1, 8)
	if _, err := r.SignalRead(nil, g1); err != nil {
		t.Fatal(err)
	}
	if _, err := r.SignalRead(nil, g2); err != nil {
		t.Fatal(err)
	}
	if !slicesEqual(g1[0][4:], g2[0][:4]) {
		t.Fatalf("overlap mismatch: g1[4:] = %v, g2[:4] = %v", g1[0][4:], g2[0][:4])
	}
}

func TestGrainReaderError(t *testing.T) {
	r := NewGrainReader(4, 2, &errReader{nch: 1})
	_, err := r.SignalRead(nil, makeGrainBuf(1, 4))
	if err == nil || err.Error() != "read error" {
		t.Fatalf("expected read error, got %v", err)
	}
}

// GrainWriter

func TestGrainWriterPrr(t *testing.T) {
	w := NewGrainWriter(4, 2, newMockWriter(1))
	sentinel := errors.New("sentinel")
	n, err := w.SignalWrite(sentinel, makeGrainBuf(1, 4))
	if n != 0 || err != sentinel {
		t.Fatalf("expected (0, sentinel), got (%d, %v)", n, err)
	}
}

func TestGrainWriterChannelPanic(t *testing.T) {
	w := NewGrainWriter(4, 2, newMockWriter(1))
	defer func() {
		if recover() == nil {
			t.Error("expected panic on channel mismatch")
		}
	}()
	w.SignalWrite(nil, makeGrainBuf(2, 4)) // writer has nch=1, grain has 2
}

func TestGrainWriterNch(t *testing.T) {
	w := NewGrainWriter(4, 2, newMockWriter(2))
	if w.NchWrite() != 2 {
		t.Fatalf("expected 2, got %d", w.NchWrite())
	}
}

func TestGrainWriterOverlapAdd(t *testing.T) {
	// nfft=4, hop=2, constant grain [1,1,1,1]
	// steady-state output per hop = grain[0]+grain[2] = 2, grain[1]+grain[3] = 2
	aw := newMockWriter(1)
	w := NewGrainWriter(4, 2, aw)
	g := makeGrainBuf(1, 4)
	for i := range g[0] {
		g[0][i] = 1
	}
	for range 4 {
		if _, err := w.SignalWrite(nil, g); err != nil {
			t.Fatal(err)
		}
	}
	// written = [0,0, 1,1, 2,2, 2,2]; steady state from index 4 onward
	want := []float64{2, 2, 2, 2}
	if !slicesEqual(aw.data[0][4:], want) {
		t.Fatalf("overlap-add: got %v, want %v", aw.data[0][4:], want)
	}
}

func TestGrainWriterError(t *testing.T) {
	w := NewGrainWriter(4, 2, &errWriter{nch: 1})
	_, err := w.SignalWrite(nil, makeGrainBuf(1, 4))
	if err == nil || err.Error() != "write error" {
		t.Fatalf("expected write error, got %v", err)
	}
}

// round-trip

func TestGrainRoundTrip(t *testing.T) {
	// hop=nfft: non-overlapping grains; writer introduces one grain of delay
	// so aw.data[0][nfft:] == input[:2*nfft] after 3 read+write cycles
	const nfft, hop = 4, 4
	input := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
	ar := &mockReader{data: [][]float64{input}, hop: hop}
	aw := newMockWriter(1)
	gr := NewGrainReader(nfft, hop, ar)
	gw := NewGrainWriter(nfft, hop, aw)

	g := makeGrainBuf(1, nfft)
	for range 3 {
		if _, err := gr.SignalRead(nil, g); err != nil {
			t.Fatal(err)
		}
		if _, err := gw.SignalWrite(nil, g); err != nil {
			t.Fatal(err)
		}
	}

	got := aw.data[0][nfft:]
	want := input[:2*nfft]
	if !slicesEqual(got, want) {
		t.Fatalf("round-trip: got %v, want %v", got, want)
	}
}
