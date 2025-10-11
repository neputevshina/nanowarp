package main

import (
	"fmt"
	"io"
	"math"
	"os"

	"image"
	"image/color"
	"image/png"

	"github.com/neputevshina/nanowarp"
	"github.com/youpy/go-wav"
)

func main() {
	nanowarp.G = map[string]any{}
	nanowarp.G[`phasogram`] = make([][]float64, 0)
	nanowarp.G[`origphase`] = make([][]float64, 0)
	nanowarp.G[`deltaphase`] = make([][]float64, 0)
	nanowarp.G[`sigmaphase`] = make([][]float64, 0)

	file, _ := os.Open(`ticktock.wav`)
	// file, _ := os.Open(`fm.wav`)
	// file, _ := os.Open(`saw-click.wav`)
	rd := wav.NewReader(file)
	data := []float64{}
	for {
		samples, err := rd.ReadSamples()
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			data = append(data, rd.FloatValue(sample, 0))
		}
	}
	nw := nanowarp.New()
	n := 2.
	out := make([]float64, int(float64(len(data)+8192)*n))
	nw.Process2(data, out, n)

	for i := range out {
		out[i] = max(-1, min(out[i], 1))
	}

	file, err := os.Create(`ticktock-x2.wav`)
	// file, err := os.Create(`fm-x2.wav`)
	// file, err := os.Create(`saw-click-x2.wav`)
	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(out)), 1, 48000, 16)
	for _, sa := range out {
		err := wr.WriteSamples([]wav.Sample{{Values: [2]int{int(sa * math.Pow(2, 16-1))}}})
		if err != nil {
			panic(err)
		}
	}
	file.Close()

	phasogram := func(name string) {
		fmt.Println(name)
		file, err := os.Create(name + `.png`)
		if err != nil {
			panic(err)
		}
		png.Encode(file, FloatMatrixToImage(nanowarp.G[name].([][]float64)))
	}
	phasogram(`phasogram`)
	phasogram(`origphase`)
	phasogram(`deltaphase`)
	phasogram(`sigmaphase`)

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
