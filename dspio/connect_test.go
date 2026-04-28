package dspio

import (
	"errors"
	"io"
	"testing"
)

func TestAltpipeNch(t *testing.T) {
	w, r := Pipe(2, 16)
	if w.NchWrite() != 2 {
		t.Errorf("NchWrite: got %d, want 2", w.NchWrite())
	}
	if r.NchRead() != 2 {
		t.Errorf("NchRead: got %d, want 2", r.NchRead())
	}
}

func TestAltpipePriorErrorShortCircuits(t *testing.T) {
	w, r := Pipe(1, 16)
	prr := errors.New("boom")
	buf := [][]float64{{1, 2, 3}}

	if n, err := w.SignalWrite(prr, buf); n != 0 || !errors.Is(err, prr) {
		t.Errorf("SignalWrite(prr): got (%d, %v), want (0, %v)", n, err, prr)
	}
	out := [][]float64{make([]float64, 3)}
	if n, err := r.SignalRead(prr, out); n != 0 || !errors.Is(err, prr) {
		t.Errorf("SignalRead(prr): got (%d, %v), want (0, %v)", n, err, prr)
	}
}

func TestAltpipeWriteChannelMismatchPanics(t *testing.T) {
	w, _ := Pipe(2, 16)
	defer func() {
		if recover() == nil {
			t.Error("expected panic on channel-count mismatch in SignalWrite")
		}
	}()
	w.SignalWrite(nil, [][]float64{{1, 2, 3}})
}

func TestAltpipeReadChannelMismatchPanics(t *testing.T) {
	_, r := Pipe(2, 16)
	defer func() {
		if recover() == nil {
			t.Error("expected panic on channel-count mismatch in SignalRead")
		}
	}()
	r.SignalRead(nil, [][]float64{make([]float64, 3)})
}

func TestAltpipeReadEmptyReturnsZero(t *testing.T) {
	_, r := Pipe(1, 16)
	out := [][]float64{make([]float64, 4)}
	n, err := r.SignalRead(nil, out)
	if err != nil {
		t.Errorf("SignalRead on empty pipe: unexpected err %v", err)
	}
	if n != 0 {
		t.Errorf("SignalRead on empty pipe: got n=%d, want 0", n)
	}
}

func TestAltpipeRoundTrip(t *testing.T) {
	w, r := Pipe(1, 16)
	in := [][]float64{{1, 2, 3, 4}}

	n, err := w.SignalWrite(nil, in)
	if err != nil {
		t.Fatalf("SignalWrite: unexpected err %v", err)
	}
	if n != len(in[0]) {
		t.Fatalf("SignalWrite: got n=%d, want %d", n, len(in[0]))
	}

	out := [][]float64{make([]float64, len(in[0]))}
	n, err = r.SignalRead(nil, out)
	if err != nil {
		t.Fatalf("SignalRead: unexpected err %v", err)
	}
	if n != len(in[0]) {
		t.Fatalf("SignalRead: got n=%d, want %d", n, len(in[0]))
	}
	for i, v := range out[0] {
		if v != in[0][i] {
			t.Errorf("out[0][%d]: got %v, want %v", i, v, in[0][i])
		}
	}
}

func TestAltpipeShortWriteWhenOverCapacity(t *testing.T) {
	w, _ := Pipe(1, 4)
	in := [][]float64{{1, 2, 3, 4, 5, 6, 7, 8}}

	n, err := w.SignalWrite(nil, in)
	if err == nil {
		t.Errorf("SignalWrite over capacity: expected error, got nil (n=%d)", n)
	}
	if err != nil && !errors.Is(err, io.ErrShortWrite) && !errors.Is(err, ErrSpillover) {
		t.Errorf("SignalWrite over capacity: got err %v, want ErrShortWrite or ErrSpillover", err)
	}
	if n > 4 {
		t.Errorf("SignalWrite over capacity: got n=%d, want <= 4", n)
	}
}

func TestAltpipeFullPipeSpills(t *testing.T) {
	w, _ := Pipe(1, 4)
	in := [][]float64{{1, 2, 3, 4}}
	if _, err := w.SignalWrite(nil, in); err != nil && !errors.Is(err, io.ErrShortWrite) {
		t.Fatalf("first write: unexpected err %v", err)
	}
	more := [][]float64{{5}}
	if _, err := w.SignalWrite(nil, more); !errors.Is(err, ErrSpillover) {
		t.Errorf("write into full pipe: got err %v, want ErrSpillover", err)
	}
}

func TestAltpipeWraparound(t *testing.T) {
	w, r := Pipe(1, 4)

	// fill, drain, then write across the wrap point
	if _, err := w.SignalWrite(nil, [][]float64{{1, 2, 3}}); err != nil {
		t.Fatalf("first write: %v", err)
	}
	out := [][]float64{make([]float64, 3)}
	if n, err := r.SignalRead(nil, out); err != nil || n != 3 {
		t.Fatalf("first read: n=%d err=%v", n, err)
	}

	// next write should wrap from index 3 to 0
	if _, err := w.SignalWrite(nil, [][]float64{{10, 20, 30}}); err != nil {
		t.Fatalf("wrap write: %v", err)
	}
	got := [][]float64{make([]float64, 3)}
	n, err := r.SignalRead(nil, got)
	if err != nil {
		t.Fatalf("wrap read: %v", err)
	}
	if n != 3 {
		t.Fatalf("wrap read: got n=%d, want 3", n)
	}
	want := []float64{10, 20, 30}
	for i, v := range got[0] {
		if v != want[i] {
			t.Errorf("wrap read got[0][%d]=%v, want %v", i, v, want[i])
		}
	}
}
