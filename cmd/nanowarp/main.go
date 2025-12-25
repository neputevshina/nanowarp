package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"runtime/pprof"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/oscope"
	"github.com/neputevshina/nanowarp/wav"
)

var println = fmt.Println

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var finput = flag.String("i", "", "input WAV (or anything else, if ffmpeg is present) `path`")
var mask = flag.Bool("mask", false, "enable auditory masking in PGHI")
var diffadv = flag.Bool("diffadv", false, "advance stereo by CIF difference, not by phase difference")
var smooth = flag.Bool("smooth", false, "trade off pre-echo for more tonal clarity")
var foutput = flag.String("o", "", "output WAV `path`")
var coeff = flag.Float64("t", 0, "time stretch multiplier")
var from = flag.Float64("from", 0, "source `bpm`")
var to = flag.Float64("to", 0, "target `bpm`")
var pitch = flag.Float64("st", 0, "pitch shift in semitones, currently adjusts time stretch without changing pitch")
var single = flag.Bool("single", false, "stretch without HPSS and using only the largest window")

func main() {
	flag.Parse()
	if *cpuprofile != "" {
		fmt.Fprintln(os.Stderr, `profiling, oscope enabled`)
		oscope.Enable = true
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

	if *coeff <= 0 && *from > 0 && *to > 0 {
		*coeff = *from / *to * math.Pow(2, *pitch/12)
	}
	if *finput == "" || *coeff <= 0 {
		flag.Usage()
		os.Exit(1)
	}
	nooutname := false
	if *foutput == "" {
		nooutname = true
		*foutput = path.Join(path.Dir(*finput), fmt.Sprintf("%.2fx-%s", *coeff, path.Base(*finput)))
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
			ex := path.Join(path.Dir(s), path.Base(s2)) + `.wav`

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

			file, err := os.Open(ex)
			if nooutname {
				*foutput = path.Join(path.Dir(ex), fmt.Sprintf("%.2fx-%s", *coeff, path.Base(ex)))
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

	mid := []float64{}
	side := []float64{}
	for {
		samples, err := wavrd.ReadSamples()
		if err == io.EOF {
			break
		}

		for _, sample := range samples {
			l, r := wavrd.FloatValue(sample, 0), wavrd.FloatValue(sample, 1)
			// mid = append(mid, l+r)
			// side = append(side, l-r)
			mid = append(mid, l)
			side = append(side, r)
		}
	}

	opts := nanowarp.Options{
		Masking: *mask,
		Smooth:  *smooth,
		Diffadv: *diffadv,
		Single:  *single,
	}
	mnw := nanowarp.New(int(wavfmt.SampleRate), opts)

	mout := make([]float64, int(float64(len(mid))**coeff))
	sout := make([]float64, int(float64(len(mid))**coeff))

	mnw.Process(mid, side, mout, sout, *coeff)

	file, err = os.Create(*foutput)

	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(mout)), 2, wavfmt.SampleRate, 32, true)
	for i := range mout {
		// msa, ssa := mout[i]/2, sout[i]/2
		// lsa, rsa := msa+ssa, msa-ssa
		lsa, rsa := mout[i], sout[i]
		err := wr.WriteSamples([]wav.Sample{{Values: [2]int{
			int(math.Float32bits(float32(lsa))),
			int(math.Float32bits(float32(rsa)))}}})
		if err != nil {
			panic(err)
		}
	}

	wd, err := os.Getwd()
	err = oscope.Dump(err, wd)
	if err != nil {
		panic(err)
	}
}
