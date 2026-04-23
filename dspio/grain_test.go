package dspio

import (
	"errors"
	"io"
	"math"
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

func TestGrainReaderOverlap(t *testing.T) {
	// nfft=8, hop=4: each consecutive grain shares nfft-hop=4 samples
	data := make([]float64, 16)
	for i := range data {
		data[i] = float64(i + 1)
	}
	r := NewGrainReader(8, 4, &mockReader{data: [][]float64{data}, hop: 4})
	g1 := make2(1, 8)
	g2 := make2(1, 8)
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

func TestGrainWriterOverlapAdd(t *testing.T) {
	// nfft=4, hop=2, constant grain [1,1,1,1]
	// steady-state output per hop = grain[0]+grain[2] = 2, grain[1]+grain[3] = 2
	aw := newMockWriter(1)
	w := NewGrainWriter(4, 2, aw)
	g := make2(1, 4)
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

	g := make2(1, nfft)
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

type phasorOsc struct {
	i, p int
}

func (p *phasorOsc) NchRead() int { return 1 }

func (p *phasorOsc) SignalRead(prr error, buf [][]float64) (n int, err error) {
	if prr != nil {
		return 0, prr
	}
	for i := range buf[0] {
		buf[0][i] = float64(p.i) / float64(p.p)
		p.i++
		p.i %= p.p
	}
	return len(buf[0]), nil
}

var _ SignalReader = &phasorOsc{}

func TestGrainAnasyn(t *testing.T) {
	// A typical analysis/synthesis scenario.
	const nfft, hop = 512, 128
	ar := &phasorOsc{p: 918}
	orig := newMockWriter(1)
	process := newMockWriter(1)
	tee := TeeReader(ar, orig)
	gr := NewGrainReader(nfft, hop, tee)
	gw := NewGrainWriter(nfft, hop, process)

	// Square of the Hann window, imitating applying it twice (before and after FFT)
	hann2 := make([]float64, nfft)
	avg := 0. // Window gain compensation
	for i := range hann2 {
		hann2[i] = 0.5 * (1 - math.Cos(2*math.Pi*float64(i)/float64(nfft)))
		hann2[i] *= hann2[i]
		avg += hann2[i]
	}
	avg /= float64(nfft)

	g := make2(1, nfft)
	for range 100 {
		n, _ := gr.SignalRead(nil, g)
		if n != hop {
			t.Fatalf("read: n != hop, got %d", n)
		}
		for i := range nfft {
			g[0][i] *= hann2[i] * avg
		}
		n, _ = gw.SignalWrite(nil, g)
		if n != hop {
			t.Fatalf("write: n != hop, got %d", n)
		}
	}

	dump("got", process.data[0], 48000)
	dump("want", orig.data[0], 48000)

	got := process.data[0][nfft:]
	want := orig.data[0][:2*nfft]

	maxerr := 0.
	for i := range min(len(got), len(want)) {
		err := got[i] - want[i]
		if math.Abs(err) > math.Abs(maxerr) {
			maxerr = err
		}
	}
	t.Log("maximum absolute error:", maxerr)

	// fmt.Println(got)

	if maxerr > 1e-4 {
		t.Fatalf("maximum error is greater than -80 dB, something is obviously wrong")
	}
}
