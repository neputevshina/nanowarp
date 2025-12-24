package main

import (
	"flag"
	"fmt"
	"math"
	"os"
)

var from = flag.Float64("from", 0, "source `bpm`")
var to = flag.Float64("to", 0, "target `bpm`")
var pitch = flag.Float64("st", 0, "pitch shift in semitones")

func main() {
	flag.Parse()
	if *from == 0 || *to == 0 {
		flag.Usage()
		os.Exit(1)
	}
	fmt.Println(*from / *to * math.Pow(2, *pitch/12))
}
