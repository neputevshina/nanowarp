package dspio

import (
	"errors"
	"fmt"
	"io"
	"sync"
	"testing"

	"github.com/neputevshina/nanowarp/waveform"
)

func TestAltpipeRoundTrip(t *testing.T) {
	w, r := lockfreePipe(1, 16)
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
	w, r := lockfreePipe(1, 4)

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

// func TestAltpipeSignal(t *testing.T) {
// 	w, r := lockfreePipe(1, 4096)

// 	wg := sync.WaitGroup{}

// 	const period = 88

// 	wg.Go(func() {
// 		b := make2(1, 2048)
// 		ar := &trainOsc{p: period}
// 		lr := LimitReader(ar, 8192*10)
// 		Copy(nil, lr, w, b)
// 		r.Close()
// 		println(`done`)
// 	})
// 	wg.Go(func() {
// 		got := make2(1, 2048)
// 		want := make2(1, 2048)
// 		etalon := &trainOsc{p: period}
// 		// for bn := 0; ; bn++ {
// 		var err error
// 		for n := 0; n != 0 || err == nil; {
// 			n, err = r.SignalRead(nil, got)
// 			if err != nil {
// 				if err == io.EOF {
// 					break
// 				}
// 				t.Fatal(err)
// 			}
// 			if n == 0 {
// 				continue
// 			}
// 			_, _ = etalon.SignalRead(nil, want)

// 			for i := range got[0] {
// 				if int(got[0][i])-int(want[0][i]) != 0 {
// 					e := waveform.Dump(nil, got[0])
// 					e = waveform.Dump(e, want[0])
// 					t.Fatal(`mismatch at`, i, `, want`, want[0][i], `got`, got[0][i])
// 				}
// 			}
// 		}
// 		// }
// 	})

// 	wg.Wait()
// }

func TestGopipeSignal(t *testing.T) {
	r, w := GoPipe(1)

	wg := sync.WaitGroup{}

	const period = 91

	wg.Add(2)
	go func() {
		b := make2(1, 2048)
		ar := &trainOsc{p: period}
		lr := LimitReader(ar, 8192*10)
		Copy(nil, lr, w, b)
		println(`done`)
		w.Close()
		wg.Done()
	}()
	erc := make(chan error)
	go func() {
		got := make2(1, 2048)
		want := make2(1, 2048)
		etalon := &trainOsc{p: period}
		// for bn := 0; ; bn++ {
		var err error
		for n := 0; n != 0 || err == nil; {
			// println(`spinout`)
			n, err = r.SignalRead(nil, got)
			if err != nil {
				if err == io.EOF {
					break
				}
				erc <- err
				return
			}
			_, _ = etalon.SignalRead(nil, want)

			for i := range got[0] {
				if int(got[0][i])-int(want[0][i]) != 0 {
					e := waveform.Dump(nil, got[0])
					e = waveform.Dump(e, want[0])
					erc <- errors.New(fmt.Sprint(`mismatch at`, i, `, want`, want[0][i], `got`, got[0][i]))
					return
				}
			}
		}
		// }
		erc <- nil
		wg.Done()
	}()
	err := <-erc
	if err != nil {
		t.Fatal(err)
	}

	wg.Wait()
}
