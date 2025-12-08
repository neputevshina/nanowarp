package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"sync"

	"image"
	"image/color"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/wav"
)

var println = fmt.Println

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var finput = flag.String("i", "", "input WAV path")
var foutput = flag.String("o", "", "output WAV path")
var coeff = flag.Float64("c", 0, "time stretch multiplier")

func main() {
	flag.Parse()
	if *cpuprofile != "" {
		fmt.Fprintln(os.Stderr, `profiling`)
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if *finput == "" || *coeff <= 0 {
		flag.Usage()
		os.Exit(1)
	}

	file, _ := os.Open(*finput)

	rd := wav.NewReader(file)
	mid := []float64{}
	side := []float64{}
	f, err := rd.Format()

	for {
		samples, err := rd.ReadSamples()
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			l, r := rd.FloatValue(sample, 0), rd.FloatValue(sample, 1)
			mid = append(mid, l+r)
			side = append(side, l-r)
		}
	}

	mnw := nanowarp.New(int(f.SampleRate))
	snw := nanowarp.New(int(f.SampleRate))

	mout := make([]float64, int(float64(len(mid))**coeff))
	sout := make([]float64, int(float64(len(mid))**coeff))

	wg := sync.WaitGroup{}
	wg.Add(2)
	go func() {
		mnw.Process(mid, mout, *coeff)
		wg.Done()
	}()
	go func() {
		snw.Process(side, sout, *coeff)
		wg.Done()
	}()
	wg.Wait()

	if *foutput == "" {
		*foutput = fmt.Sprintf("%.2fx-%s", *coeff, *finput)
	}

	file, err = os.Create(*foutput)

	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(mout)), 2, f.SampleRate, 32, true)
	for i := range mout {
		msa, ssa := mout[i]/2, sout[i]/2
		lsa, rsa := msa+ssa, msa-ssa
		err := wr.WriteSamples([]wav.Sample{{Values: [2]int{
			int(math.Float32bits(float32(lsa))),
			int(math.Float32bits(float32(rsa)))}}})
		if err != nil {
			panic(err)
		}
	}
	file.Close()
}

// FloatMatrixToImage converts a 2D float64 matrix into a grayscale image.Image.
// The function normalizes the input values to the range [0, 255].
//
// Courtesy of ChatGPT.
func FloatMatrixToImage(data [][]float64) image.Image {
	if len(data) == 0 || len(data[0]) == 0 {
		return nil
	}

	height := len(data[0])
	width := len(data)

	// Find min and max
	minVal := math.Inf(1)
	maxVal := math.Inf(-1)
	for _, row := range data {
		for _, v := range row {
			if v < minVal {
				minVal = v
			}
			if v > maxVal {
				maxVal = v
			}
		}
	}
	fmt.Println(minVal, maxVal)

	scale := 1.
	offset := 3.14
	scale = 255.0 / (maxVal - minVal)
	offset = -minVal * scale

	img := image.NewGray(image.Rect(0, 0, width, height))
	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			val := data[x][y]*scale + offset
			if val < 0 {
				val = 0
			}
			if val > 255 {
				val = 255
			}
			img.SetGray(x, y, color.Gray{Y: uint8(val + 0.5)})
		}
	}

	return img
}
