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
	"os/signal"
	"path"
	"path/filepath"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"slices"
	"strconv"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/dspio/wavio"
	"github.com/neputevshina/nanowarp/oscope"
)

var println = fmt.Println

var progress = flag.Bool("p", true, "Display progress bar.")
var cpuprofile = flag.String("cpuprofile", "", "Write cpu profile to `file`.")
var finput = flag.String("i", "", "Input WAV (or anything else, if ffmpeg is present) `path`.")
var foutput = flag.String("o", "", "Output WAV `path`.")
var from = flag.Float64("from", 1, "Source `BPM`.")
var to = flag.Float64("to", 1, "Target `BPM` or stretch factor if -from is not set.")
var st = flag.Float64("st", 0, `Pitch shift in semitones.
Currently adjusts time stretch without changing the the pitch.`)
var onsets = flag.Bool("onsets", false, "Output displaced onsets only.")
var resets = flag.Int("resets", 0, `Time and phase resets:
-2: Don't perform transient separation, output raw PVDR without phase resets.
-1: Use original phase when not stretching.
    Introduces clicky artifacts but cleanest for transient-heavy material.
    Best numerical stability because of resets.
0:  When not stretching, advance phase for harmonic components, 
    use original phase otherwise.
    Very little to no artifacts, but noticeable slight loss in clarity.`)
var q = flag.Int("q", 0, `Quality:
Set algorithm quality.
-1: Use brute force approximation to PGHI. Less transparent, 20% faster.
0:  Use PGHI.`)
var poolms = flag.Int("poolms", 200, `Time of onset detection bucket in milliseconds.
Minimum amount of time between two consecutive transient detections.`)
var outpool = flag.Bool("outpool", false, `If true, measure pooling size in output time, not in input time.
I.e. scale the pooling size with the stretch coefficient.
Effective only for stretches, not shrinks, which are always scaled.`)
var ri = flag.Int("ri", 10, `Minimum amount of time bins a trajectory must travel in total
to be considered tonal.
Lower values — more tonal preservation and less transient clarity.
Higher values — more transient preservation and more interrupts.`)
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
var notriple = flag.Int("triplefix", 1, `Behavior of limiting of local group delay.`)
var nfft = flag.Int("nfft", 3600, "Base window size at 48000 Hz sample rate.")

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

Basic usage:
  -i path
	Input WAV (or anything else, if ffmpeg is present) path.
  -o path
	Output WAV path. If empty will be generated from input path.
  -from BPM
	Source BPM. (default 1)
  -to BPM
	Target BPM or stretch factor if -from is 1. (default 1)
  -st float
	Pitch shift in semitones.
	Currently adjusts time stretch without changing the the pitch.
  -timemap string
	Path to timemap file in rubberband program format.
	Each line is a pair of two integers, separated by space.
	First integer is an input sample index, second is an output sample index.
	Unlike rubberband, specifying total duration is not needed, but last pair
	must be a pair of input sample index and a last output sample index.
	Time map must be functional. Output index of any breakpoint can't be
	less than that of any previous breakpoint.

Quality and behavior:
  -q int
	Set algorithm quality.
	-1: Use brute force approximation to PGHI. Less transparent, 20% faster.
	0:  Use PGHI.
  -resets int
	Time and phase resets:
	-2: Don't perform transient separation, output raw PVDR without phase resets.
	-1: Use original phase when not stretching.
	    Introduces clicky artifacts but cleanest for transient-heavy material.
	    Best numerical stability because of resets.
	0:  When not stretching, advance phase for harmonic components,
	    use original phase otherwise.
	    Very little to no artifacts, but noticeable slight loss in clarity.
  -nfft int
	Base window size in samples at 48000 Hz. Discriminates between bass and pulsation.
	Lower values give more clarity at the cost of increasing lowest possible note.
	The correct value tends to correspond to half of the frequency
	of the lowest fundamental in the signal.
	Try values in range 3000–4096.
	Values greater than 4096 will oversample the FFT 4 or more times,
	increasing run time __without__ increase in quality. 
	Default is 3600.
  -outpool
	If true, measure pooling size in output time, not in input time.
	I.e. scale the pooling size with the stretch coefficient.
	Effective only for stretches, not shrinks, which are always scaled.  
  -poolms int
	Time of onset detection bucket in milliseconds.
	Minimum amount of time between two consecutive transient detections. 
	Default is 200.
  -ri int
	Minimum amount of time bins a trajectory must travel in total
	to be considered tonal.
	Lower values — more tonal preservation and less transient clarity.
	Higher values — more transient preservation and more interrupts.
	Default is 10.
  -if int
	Maximum radius of influence of each detected tonal trajectory.
	Phase never be reset at this number of bins around the ridge.
	Higher values compromise transient quality over tonal quality.
	Default is 2.
  -triplefix
	Behavior of limiting of local group delay.
	-1 disables fix for ridge triplication in time, which is obvious on
	extreme (>4x) stretches. 
	0 enables it for coefficients strictly greater than 2.
	1 forces it. 
	Default is 1.

Utility:
  -p    Display progress bar. (default true)

Debug:
  -cpuprofile file
	Write cpu profile to file.
  -onsets
	Output displaced onsets only.
  -experiment number
	DON'T USE: run a numbered experiment instead of nanowarp.
  -flac string
	Not implemented. Output FLAC encoded file.
`)
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
	generateOutSuffix := func() string {
		pitchSuffix := ""
		if *st != 0 {
			pitchSuffix = fmt.Sprintf("%+.2fst", *st)
		}

		if *from > 1 {
			return fmt.Sprintf("%g→%g%s", *from, *to, pitchSuffix)
		} else if math.Abs(*st) > 0 {
			return fmt.Sprintf("%s", pitchSuffix)
		} else {
			return fmt.Sprintf("%.4fx%s", *to, pitchSuffix)
		}
	}
	generateOutName := func(dir, fn string) string {
		return path.Join(path.Dir(dir), fmt.Sprintf("%s-%s", generateOutSuffix(), path.Base(fn)))
	}
	if *foutput == "" {
		nooutname = true
		*foutput = generateOutName(path.Dir(*finput), path.Base(*finput))
	}

	file, err := os.Open(*finput)
	if err != nil {
		panic(err)
	}

	wsr, err := wavio.NewDecoder(file)
	if err != nil {
		// Try to call ffmpeg, we've probably got an MP3.
		if err == wavio.ErrNotAWav {
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

			cmd := exec.Command(`ffmpeg`, `-hide_banner`, `-y`, `-i`, s2, `-acodec`, `pcm_f32le`, ex)
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
			wsr, err = wavio.NewDecoder(file)
			if err != nil {
				panic(err)
			}
		} else {
			panic(err)
		}
	}

	props := wsr.Properties()

	inputLength := float64(props.Samples)

	if *experiment != 0 {
		// experiments(*experiment, file, *foutput)
		panic(`no experiments today`)
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
		if *from > 1 || math.Abs(*st) > 0 {
			*to = *from / *to * math.Pow(2, *st/12)
		}
		if *to <= 0 {
			flag.Usage()
			os.Exit(1)
		}
		bps = []nanowarp.Breakpoint{nanowarp.Bp(0, 0), nanowarp.Bp(inputLength, inputLength**to)}
	}
	slices.SortFunc(bps, func(a, b nanowarp.Breakpoint) int {
		return int((a.I - b.I) / math.Abs(a.I-b.I))
	})
	phasor, err := nanowarp.NewCurve(bps)
	if err != nil {
		panic(err)
	}

	// Open the temporary output and delete it if interrupted.
	outfile, err := os.CreateTemp(path.Dir(*foutput), "nanowarp_*.wav")
	untmp := func() {
		_ = outfile.Close()
		_ = os.Remove(outfile.Name())
	}
	defer untmp()
	if err != nil {
		panic(err)
	}
	irq := make(chan os.Signal, 1)
	signal.Notify(irq, os.Interrupt)
	go func() {
		<-irq
		untmp()
		os.Exit(0)
	}()

	wsw, err := wavio.NewEncoder(outfile, props.Samplerate, props.Nch, wavio.FormatFloat, 32)
	if err != nil {
		panic(err)
	}

	// Copy all sections from original to the output.
	for _, hd := range wsr.Headermap {
		if string(hd.FourCC[:]) == `fmt ` || string(hd.FourCC[:]) == `data` {
			continue
		}
		_, err := file.Seek(hd.Seek, io.SeekStart)
		if err != nil {
			panic(err)
		}
		err = wsw.WriteRiffChunk(hd.FourCC, io.LimitReader(file, int64(hd.Size)))
		if err != nil {
			panic(err)
		}
	}

	wsr.Rewind()

	// Warping.
	var pch chan nanowarp.Progress
	if *progress {
		pch = make(chan nanowarp.Progress)
	}
	opts := nanowarp.Options{
		Onsets:   *onsets,
		Quality:  *q,
		Resets:   *resets,
		Progress: pch,
		Hyperparams: nanowarp.Hyperparams{
			PickingMs:       *poolms,
			ScalePool:       *outpool,
			InfluenceRadius: *ifr,
			LongRidgeLength: *ri,
			TriplicationFix: *notriple,
			BaseNFFT:        *nfft,
		},
	}
	tsm := nanowarp.New(props.Samplerate, props.Nch, opts)

	var exit chan struct{}
	if *progress {
		exit = make(chan struct{})
		pb := startProgress(os.Stderr)
		go func() {
			for bp := range pch {
				pb.Set(bp.Current, bp.End)
				fmt.Fprint(os.Stderr, " ", bp.Process)
			}
			fmt.Fprintln(os.Stderr)
			exit <- struct{}{}
		}()
	}

	tsm.Process(&withlength{wsr}, wsw, phasor)

	err = wsw.Close()
	if err != nil {
		panic(err)
	}
	err = file.Close()
	if err != nil {
		panic(err)
	}

	if exit != nil {
		<-exit
	}

	// Moving the temporary output to real one.
	err = os.Rename(outfile.Name(), *foutput)
	if err != nil {
		panic(err)
	}
	oscope.Dump(nil, "./pics")
}

type withlength struct {
	*wavio.Decoder
}

func (w *withlength) Length() int { return w.Properties().Samples }
