package main

import (
	"fmt"
	"io"
	"math"
	"os"
	"sync"

	"image"
	"image/color"

	"github.com/neputevshina/nanowarp"
	"github.com/youpy/go-wav"
)

var println = fmt.Println

func main() {
	// filename := `fm.wav`
	// filename := `output.wav`
	// filename := `saw.wav`
	// filename := `saw-short.wav`
	// filename := `saw-click.wav`
	// filename := `saw-click.wav`
	// filename := `ticktock.wav`
	// filename := `welcome.wav`
	filename := `ЗАПАХЛО_NIGHTCALL_МЭШАПЕР_АРКАДИЙ_ГАЧИБАСОВ.mp3.wav`
	// filename := `Диалоги тет-а-тет - ALEKS ATAMAN.m4a.mp3.wav`
	// filename := `audio_2025-12-04_04-07-32.ogg.wav`

	fmt.Fprintln(os.Stderr, filename)

	file, _ := os.Open(filename)

	rd := wav.NewReader(file)
	left := []float64{}
	right := []float64{}
	f, err := rd.Format()

	for {
		samples, err := rd.ReadSamples()
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			left = append(left, rd.FloatValue(sample, 0))
			right = append(right, rd.FloatValue(sample, 1))
		}
	}

	lnw := nanowarp.New(int(f.SampleRate))
	rnw := nanowarp.New(int(f.SampleRate))

	var n float64 = 2
	lout := make([]float64, int(float64(len(left)+8192)*n))
	rout := make([]float64, int(float64(len(left)+8192)*n))

	wg := sync.WaitGroup{}
	wg.Add(2)
	go func() {
		lnw.Process(left, lout, n)
		wg.Done()
	}()
	go func() {
		rnw.Process(right, rout, n)
		wg.Done()
	}()
	wg.Wait()

	for i := range lout {
		lout[i] *= 0.25
		rout[i] *= 0.25
	}

	file, err = os.Create(fmt.Sprintf("%.2fx-%s", n, filename))

	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(lout)), 2, f.SampleRate, 16)
	for i := range lout {
		lsa, rsa := lout[i], rout[i]
		err := wr.WriteSamples([]wav.Sample{{Values: [2]int{int(lsa * math.Pow(2, 16-1)), int(rsa * math.Pow(2, 16-1))}}})
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
