package dspio

import (
	"io"
	"sync"
	"testing"

	"github.com/neputevshina/nanowarp/waveform"
)

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

func TestAltpipeSignal(t *testing.T) {
	w, r := Pipe(1, 4096)

	wg := sync.WaitGroup{}

	wg.Go(func() {
		b := make2(1, 2048)
		ar := &trainOsc{p: 55}
		lr := LimitReader(ar, 8192*10)
		Copy(nil, lr, w, b)
	})
	wg.Go(func() {
		got := make2(1, 2048)
		want := make2(1, 2048)
		etalon := &trainOsc{p: 55}
		for bn := 0; ; bn++ {
			_, err := r.SignalRead(nil, got)
			if err == io.EOF {
				break
			}
			if err != nil {
				t.Fatal(err)
			}
			_, _ = etalon.SignalRead(nil, want)
			for i := range got[0] {
				if int(got[0][i])-int(want[0][i]) != 0 {
					waveform.Dump(nil, got[0])
					waveform.Dump(nil, want[0])
					t.Fatal(`mismatch at`, i, `, bn`, bn, `want`, want[0][i], `got`, got[0][i])
				}
			}
		}
	})

	wg.Wait()
}
