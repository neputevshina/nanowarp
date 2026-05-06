// Package waveform provides a way to dump the waveform view of a buffer to the terminal.
package waveform

import (
	"fmt"
	"os"
	"slices"

	"golang.org/x/term"
)

func Dump(prr error, b []float64) error {
	const height = 7
	if prr != nil {
		return prr
	}
	w, _, err := term.GetSize(1)
	if err != nil {
		return err
	}
	bin := max(1, len(b)/w)
	gmax, gmin := slices.Max(b), -slices.Min(b)
	absmax := max(gmax, -gmin)
	if absmax != 0 {
		for i := range b {
			b[i] /= absmax
		}
	}

	lines := make([][]byte, height)
	column := make([]byte, height)
	for i := 0; i < len(b); i += bin {
		fill(column, ' ')
		sl := b[i:min(len(b), i+bin)]
		ma, mi := (-slices.Max(sl))/2+0.5, (-slices.Min(sl))/2+0.5
		h := float64(height)
		column[int(h*mi)-1] = '\''
		// fill(column[int(h*ma):int(h*mi)], ':')
		column[int(h*ma)] = '.'
		for i := range height {
			lines[i] = append(lines[i], column[i])
		}
	}

	sa := fmt.Sprint(gmax)
	f := fmt.Sprintf("%%s%%%dd\n", w-len(sa))
	fmt.Fprintf(os.Stderr, f, sa, len(b))
	for _, l := range lines {
		_, err := os.Stderr.Write(l[:w])
		if err != nil {
			return err
		}
		_, err = os.Stderr.Write([]byte{'\n'})
		if err != nil {
			return err
		}
	}
	fmt.Fprintln(os.Stderr, gmin)

	return nil
}

func fill[T any](s []T, e T) {
	for i := range s {
		s[i] = e
	}
}
