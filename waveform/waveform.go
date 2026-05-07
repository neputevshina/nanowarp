// Package waveform provides a way to dump the waveform view of a buffer to the terminal.
package waveform

import (
	"fmt"
	"math"
	"os"
	"slices"

	"golang.org/x/term"
)

func Dump(prr error, b []float64) error {
	const height = 12
	if prr != nil {
		return prr
	}
	w, _, err := term.GetSize(1)
	if err != nil {
		return err
	}
	symb := max(1, len(b)/w)
	gmax, gmin := slices.Max(b), slices.Min(b)
	absmax := max(gmax, -gmin)
	if absmax != 0 {
		for i := range b {
			b[i] /= absmax
		}
	}

	lines := make([][]rune, height)
	column := make([]rune, height)
	for i := 0; i < len(b); i += symb {
		fill(column, ' ')
		sl := b[i:min(len(b), i+symb)]
		h := float64(height)

		floor := func(x float64) int { return int(math.Round(x)) }
		frac := func(x float64) float64 { _, f := math.Modf(x); return f }

		ma, mi := h*((-slices.Max(sl))/2+0.5), h*((-slices.Min(sl))/2+0.5)
		if ma != mi {
			fill(column[floor(ma):min(height, floor(mi)+1)], '█')
		}
		if frac(ma) <= 0.5 {
			column[min(height-1, floor(mi))] = '▀'
		}
		if frac(mi) >= 0.5 {
			column[floor(ma)] = '▄'
		}

		for i := range height {
			lines[i] = append(lines[i], column[i])
		}
	}

	sa := fmt.Sprint(gmax)
	f := fmt.Sprintf("%%s%%%dd\n", w-len(sa))
	fmt.Fprintf(os.Stderr, f, sa, len(b))
	for _, l := range lines {
		_, err := os.Stderr.Write([]byte(string(l[:w])))
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
