package main

import (
	"fmt"
	"math"
	"os"
	"time"

	"golang.org/x/term"
)

type Progress struct {
	f            *os.File
	Prev         float64
	Start, Prevt time.Time
}

var semi = []rune(" ▏▎▍▌▋▊▉█")

func startProgress(w *os.File) *Progress {
	return &Progress{f: w, Start: time.Now()}
}

func (p *Progress) Set(nu, last float64) {
	n := nu / last
	w, _, _ := term.GetSize(1) // Don't care about errors
	ew := w / 2
	fmt.Fprint(p.f, "\r")
	for range w {
		fmt.Fprint(p.f, " ")
	}
	i, f := math.Modf(n * float64(ew))
	fmt.Fprintf(p.f, "\r% 4.f%% ", n*100)
	for range int(i) {
		fmt.Fprint(p.f, string(semi[len(semi)-1]))
	}
	fmt.Fprint(p.f, string(semi[int(f*float64(len(semi)))]))
	for range ew - int(i) {
		fmt.Fprint(p.f, string(semi[0]))
	}

	now := time.Now()
	elapsed := now.Sub(p.Start).Round(time.Second)
	// eta := float64(elapsed) / n * p.Last
	fmt.Fprint(p.f, elapsed)

	p.Prev = nu
	p.Prevt = now
}
