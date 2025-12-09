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
	"sync"

	"github.com/neputevshina/nanowarp"
	"github.com/neputevshina/nanowarp/wav"
)

var println = fmt.Println

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var finput = flag.String("i", "", "input WAV (or anything else, if ffmpeg is present) `path`")
var foutput = flag.String("o", "", "output WAV `path`")
var coeff = flag.Float64("t", 0, "time stretch multiplier")

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
			mid = append(mid, l+r)
			side = append(side, l-r)
		}
	}

	mnw := nanowarp.New(int(wavfmt.SampleRate))
	snw := nanowarp.New(int(wavfmt.SampleRate))

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

	file, err = os.Create(*foutput)

	if err != nil {
		panic(err)
	}
	wr := wav.NewWriter(file, uint32(len(mout)), 2, wavfmt.SampleRate, 32, true)
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
