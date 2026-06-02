// Package waveform provides a way to dump the pseudographic waveform view of an audio buffer.
package waveform

import (
	"fmt"
	"io"
	"math"
	"os"
	"slices"

	"golang.org/x/exp/constraints"
	"golang.org/x/term"
)

// Lookup is the lookup table used to draw a 2×3 dot rectangle with symbols.
//
// The order of bits in the rectangle of dots by which it is indexed is:
//
//	0 1
//	2 3
//	4 5
//
// It must contain at least 64 elements.
var Lookup = []rune(` 🬀🬁🬂🬃🬄🬅🬆🬇🬈🬉🬊🬋🬌🬍🬎🬏🬐🬑🬒🬓▌🬔🬕🬖🬗🬘🬙🬚🬛🬜🬝🬞🬟🬠🬡🬢🬣🬤🬥🬦🬧▐🬨🬩🬪🬫🬬🬭🬮🬯🬰🬱🬲🬳🬴🬵🬶🬷🬸🬹🬺🬻█`) // U+1FB00-U+1FB3B is a lie.

// Dump prints the waveform of the audio to Stderr.
func Dump(prr error, b []float64) error {
	return DumpWriter(prr, os.Stderr, b)
}

// DumpWrited prints the waveform of the audio to out.
func DumpWriter(prr error, out io.Writer, b []float64) error {
	b = slices.Clone(b)
	const height = 12
	const pixheight = height * 3
	if prr != nil {
		return prr
	}
	w, _, err := term.GetSize(1)
	if err != nil {
		return err
	}
	w = min(w, len(b))
	symb := max(1, len(b)/w/2)
	gmax, gmin := slices.Max(b), slices.Min(b)
	absmax := max(gmax, -gmin)
	if absmax != 0 {
		for i := range b {
			b[i] = unmix(gmin, gmax, b[i])
		}
	}

	column := make([]bool, pixheight)
	blines := make([][]bool, pixheight)
	for i := 0; i < len(b); i += symb {
		fill(column, false)
		sl := b[i:min(len(b), i+symb)]
		h := float64(pixheight)

		floor := func(x float64) int { return int(math.Floor(x)) }
		ceil := func(x float64) int { return int(math.Ceil(x)) }

		ma := clamp(0, len(column)-1, floor(h*(1-slices.Max(sl))))
		mi := clamp(0, len(column)-1, ceil(h*(1-slices.Min(sl))))
		fill(column[ma:mi], true)
		if ma == mi {
			column[ma] = true
		}

		for i := range pixheight {
			blines[i] = append(blines[i], column[i])
		}
	}
	g := func(x, y int) int {
		if y >= len(blines) {
			return 0
		}
		if x >= len(blines[0]) {
			return 0
		}
		if blines[y][x] {
			return 1
		}
		return 0
	}
	lines := make([][]rune, height)
	for i := range lines {
		lines[i] = make([]rune, w)
	}
	for x := 0; x < len(blines[0]); x += 2 {
		for y := 0; y < len(blines); y += 3 {
			i := g(x, y) | g(x+1, y)<<1 |
				g(x, y+1)<<2 | g(x+1, y+1)<<3 |
				g(x, y+2)<<4 | g(x+1, y+2)<<5
			r := Lookup[i]
			if y/3 < len(lines) && x/2 < len(lines[0]) {
				lines[y/3][x/2] = r
			}
		}
	}

	sa := fmt.Sprint(gmax)
	f := fmt.Sprintf("%%s%%%dd\n", w-len(sa))
	fmt.Fprintf(out, f, sa, len(b))
	for _, l := range lines {
		_, err := out.Write([]byte(string(l[:w])))
		if err != nil {
			return err
		}
		_, err = out.Write([]byte{'\n'})
		if err != nil {
			return err
		}
	}
	fmt.Fprintln(out, gmin)

	return nil
}

func fill[T any](s []T, e T) int {
	for i := range s {
		s[i] = e
	}
	return len(s)
}

func unmix[F constraints.Float](a, b, x F) F {
	return (x - a) / (b - a)
}

func clamp[T constraints.Ordered](a, b, x T) T {
	return max(a, min(b, x))
}
