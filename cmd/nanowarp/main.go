package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"slices"
	"strconv"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/wav"
	"github.com/schollz/progressbar/v3"
)

var println = fmt.Println

var progress = flag.Bool("p", false, "Display progress bar.")
var cpuprofile = flag.String("cpuprofile", "", "Write cpu profile to `file`.")
var finput = flag.String("i", "", "Input WAV (or anything else, if ffmpeg is present) `path`.")
var foutput = flag.String("o", "", "Output WAV `path`.")
var coeff = flag.Float64("t", 0, "Time stretch multiplier.")
var from = flag.Float64("from", 1, "Source `BPM`.")
var to = flag.Float64("to", 1, "Target `BPM`.")
var st = flag.Float64("st", 0, `Pitch shift in semitones.
Currently adjusts time stretch without changing the the pitch.`)
var onsets = flag.Bool("onsets", false, "Output displaced onsets only.")
var q = flag.Int("q", 0, `Quality:
-2: Don't perform transient separation, output raw PVDR without phase resets.
-1: Extract transients and reset the phase when not stretching.
    Introduces clicky artifacts but cleanest for transient-heavy material.
    Best numerical stability because of resets.
0:  Same as -1, but detects and bypasses tonal components.
    No artifacts, but noticeable slight loss in clarity.`)
var poolms = flag.Int("poolms", 250, `Time of onset detection bucket in milliseconds.
Minimum amount of time between two consecutive transient detections.`)
var outpool = flag.Bool("outpool", false, `If true, measure pooling size in output time, not in input time.
I.e. scale the pooling size with the stretch coefficient.
Effective only for stretches, not shrinks, which are always scaled.`)
var ifr = flag.Int("if", 2, `Maximum radius of influence of each detected tonal trajectory.
Phase never be reset at this number of bins around the ridge.
Higher values compromise transient quality over tonal quality.`)
var timemappath = flag.String("timemap", "", `Path to timemap file in rubberband program format.
Each line is a pair of two integers, separated by space. 
First integer is an input sample index, second is an output sample index.
Unlike rubberband, specifying total duration is not needed, but last pair 
must be a pair of input sample index and a last output sample index.
Time map must be functional. Output index of any breakpoint can't be 
less than that of any previous breakpoint.`)
var flac = flag.String("flac", "", `Not implemented. Output FLAC encoded file.`)
var experiment = flag.Int("experiment", 0, "DON'T USE: run a `number`ed experiment instead of nanowarp.")

func init() {
	flag.Usage = func() {
		buildinfo, _ := debug.ReadBuildInfo()
		fmt.Fprint(flag.CommandLine.Output(),
			`Nanowarp `,
			buildinfo.Main.Version, ` `, runtime.Version(),
			`
Audio time stretching algorithm
© 2025-2026 neputevshina
https://github.com/neputevshina/nanowarp

Usage:
`)
		flag.PrintDefaults()
	}
	flag.Parse()
}

func main() {
	var err error
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

	if *finput == "" {
		flag.Usage()
		os.Exit(1)
	}

	nooutname := false
	generateOutName := func(dir, fn string) string {
		pitchSuffix := ""
		if *st != 0 {
			pitchSuffix = fmt.Sprintf("%+.2fst", *st)
		}

		if *from > 1 {
			return path.Join(path.Dir(dir), fmt.Sprintf("%g→%g%s-%s", *from, *to, pitchSuffix, path.Base(fn)))
		} else if math.Abs(*st) > 0 {
			return path.Join(path.Dir(dir), fmt.Sprintf("%s-%s", pitchSuffix, path.Base(fn)))
		} else {
			return path.Join(path.Dir(dir), fmt.Sprintf("%.4fx%s-%s", *coeff, pitchSuffix, path.Base(fn)))
		}
	}
	if *foutput == "" {
		nooutname = true
		*foutput = generateOutName(path.Dir(*finput), path.Base(*finput))
	}

	file, err := os.Open(*finput)
	if err != nil {
		panic(err)
	}

	wavrd := wav.NewReader(file)
	wavfmt, err := wavrd.Format()
	if err != nil {
		// Try to call ffmpeg, we've probably got an MP3.
		// FIXME Bodge error type check
		if err.Error() == `Given bytes is not a RIFF format` {
			_ = file.Close()

			_, err = exec.LookPath(`ffmpeg`)
			if err != nil {
				fmt.Fprintln(os.Stderr, `can process only WAV files without ffmpeg`)
				os.Exit(1)
			}

			s, err := filepath.Abs(*foutput)
			if err != nil {
				panic(err)
			}
			s2, err := filepath.Abs(*finput)
			if err != nil {
				panic(err)
			}

			// ex := path.Join(path.Dir(s), path.Base(s2)) + `.wav`
			ex := path.Join(os.TempDir(), path.Base(s2)) + `.wav`

			cmd := exec.Command(`ffmpeg`, `-y`, `-i`, s2, `-acodec`, `pcm_f32le`, ex)
			cmd.Stderr = os.Stderr
			cmd.Stdout = os.Stdout
			err = cmd.Run()
			if err != nil {
				if _, ok := err.(*exec.ExitError); ok {
					os.Exit(1)
				} else {
					panic(err)
				}
			}

			file, err = os.Open(ex)
			if nooutname {
				*foutput = generateOutName(s, ex)
			}
			if err != nil {
				panic(err)
			}
			wavrd = wav.NewReader(file)
			wavfmt, err = wavrd.Format()
			if err != nil {
				panic(err)
			}
		} else {
			panic(err)
		}
	}

	inputLength := 0.
	_, err = wavrd.Duration() // Populate wavrd.Size field
	if err != nil {
		panic(err)
	}
	inputLength = float64(wavrd.Size) / float64(wavfmt.BlockAlign)
	// TODO This wav package is a piece of shit. Write own.

	if *experiment != 0 {
		experiments(*experiment, file, *foutput)
		return
	}

	// Phase ramp generation.
	var bps []nanowarp.Breakpoint
	var timemapfile *os.File
	if *timemappath != "" {
		timemapfile, err = os.Open(*timemappath)
		if err != nil {
			panic(err)
		}
		bs, err := io.ReadAll(timemapfile)
		if err != nil {
			panic(err)
		}
		lines := bytes.SplitSeq(bs, []byte("\n"))
		for l := range lines {
			n := len(bps)
			ib, jb, ok := bytes.Cut(l, []byte(` `))
			if !ok {
				panic(fmt.Sprint(`malformed index pair at line `, n))
			}
			// 53 is the lossless integer resolution of float64 in bits.
			// In an unlikely event some of indexes got greater
			// than 9 quadrillion (which at 48000 Hz is ≈6000 years),
			// strconv will catch it.
			i, err := strconv.ParseInt(string(ib), 10, 53)
			if err != nil {
				panic(fmt.Sprint(`first element of pair is invalid number, line `, n))
			}
			j, err := strconv.ParseInt(string(jb), 10, 53)
			if err != nil {
				panic(fmt.Sprint(`second element of pair is invalid number, line `, n))
			}
			bps = append(bps, nanowarp.Bp(float64(i), float64(j)))
		}
	} else {
		if *coeff <= 0 && (*from > 1 && *to > 1 || math.Abs(*st) > 0) {
			*coeff = *from / *to * math.Pow(2, *st/12)
		}
		if *coeff <= 0 {
			flag.Usage()
			os.Exit(1)
		}
		bps = []nanowarp.Breakpoint{nanowarp.Bp(0, 0), nanowarp.Bp(inputLength, inputLength**coeff)}
	}
	slices.SortFunc(bps, func(a, b nanowarp.Breakpoint) int {
		return int((a.I - b.I) / math.Abs(a.I-b.I))
	})
	phasor, err := nanowarp.NewCurve(bps)
	if err != nil {
		panic(err)
	}

	// Reading the whole file to memory.
	left := []float64{}
	right := []float64{}
	var samples []wav.Sample
	for {
		_, samples, err = wavrd.ReadSamples(samples)
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			l, r := wavrd.FloatValue(sample, 0), wavrd.FloatValue(sample, 1)
			left = append(left, l)
			right = append(right, r)
		}
	}

	// Working.
	var pch chan nanowarp.Breakpoint
	if *progress {
		pch = make(chan nanowarp.Breakpoint)
	}
	opts := nanowarp.Options{
		Onsets:   *onsets,
		Quality:  *q,
		Progress: pch,
		Hyperparams: nanowarp.Hyperparams{
			PickingMs:       *poolms,
			ScalePool:       *outpool,
			InfluenceRadius: *ifr,
		},
	}
	mnw := nanowarp.New(int(wavfmt.SampleRate), opts)

	end := int(bps[len(bps)-1].J)
	lout := make([]float64, end)
	rout := make([]float64, end)

	if *progress {
		pb := progressbar.NewOptions(end,
			progressbar.OptionShowElapsedTimeOnFinish(),
			progressbar.OptionSetTheme(progressbar.Theme{
				Saucer:        "█",
				SaucerPadding: " ",
			}))
		go func() {
			for bp := range pch {
				pb.Set(int(bp.J))
			}
			println()
		}()
	}
	mnw.Process(left, right, lout, rout, phasor)

	// Dumping the output.
	file, err = os.Create(*foutput)
	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(lout)), 2, wavfmt.SampleRate, 32, true)

	nbuf := 2048
	buf := make([]wav.Sample, 0, nbuf)
	for i := 0; i < len(lout); i += nbuf {
		buf = buf[:0]
		for j := i; j < min(i+nbuf, len(lout)); j++ {
			lsa, rsa := lout[j], rout[j]
			buf = append(buf, wav.Sample{Values: [2]int{
				int(math.Float32bits(float32(lsa))),
				int(math.Float32bits(float32(rsa)))}})
		}
		err := wr.WriteSamples(buf)
		if err != nil {
			panic(err)
		}
	}

	oscope.Dump(nil, "./pics")
}
