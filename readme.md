# Nanowarp
Studio-grade audio time stretching algorithm.

Root package is Nanowarp library, cmd/nanowarp is CLI.

## Installation

```
go get github.com/neputevshina/nanowarp
```

## nanowarp CLI: Installation and usage

1. Install [Go](https://go.dev/)
2. Open terminal/shell/command prompt.
3. Install Nanowarp
```
go install github.com/neputevshina/nanowarp/cmd/nanowarp@latest
```
Please note that `@latest` installs the last tagged version (e.g. v0.5.4), not the latest commit in master.

4. Use it
```
nanowarp -i inputfile.wav -to <stretch> [-o outputfile.wav]
```
or
```
nanowarp -i inputfile.wav -from <bpm> -to <bpm> -st <semitones> [-o outputfile.wav]
```
or
```
nanowarp -i inputfile.wav -timemap <time map file> [-o outputfile.wav]
```
If your system can't find `nanowarp` executable, you have probably changed PATH environment variable.
Probably the simplest way to bring it back if you are under Windows is by reinstalling the Go.
On Linux and macOS, you should know what to do.

Consult 
```
nanowarp -help
```
to get the list of available options.

## Implementation

Nanowarp is a phase gradient heap integration (PGHI) phase vocoder (aka PVDR)[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two for simplicity.

A phase ramp for the entire output signal is generated. Onsets are detected using a variant 
of Superflux method[3]. If onset is detected, phase ramp will have a derivative (scan speed) of 
1 in a region around detected onset. Starting points of these sample regions are scaled by the 
stretch size, and points between regions are linearly interpolated.

Then the large-grained (nfft=4096) PVDR is applied, using phase ramp curve for input 
sample indexes. PVDR is applied in two steps: first, integration directions (arrows) are detected.
If `Quality` parameter set to -1, Nanowarp uses faster greedy local maximum arrow detection instead of PGHI, which uses priority queue.
Sources of these arrows (ridges) are extracted too. Then, phase is accumulated from partial derivatives by interpreting the arrows. 
If current frame speed is 1, phase accumulator is reset, except ridges from previous 
frame with their regions of influence, which are accumulated further.

Like in original implementation of PVDR, FFT is oversampled by factor of 2 with zero-padding. 
Stereo coherence is obtained through stretching mono and adding complex phase difference of 
respective side channels after stretching back[4].

Transient triplication, characteristic to all PVDRs, is mitigated by removing all bins where 
absolute horizontal displacement (partial derivative of phase by time) is greater than 
half analysis hop size. This mechanism is engaged only for stretch coefficients greater than 2.

The algorithm does not depend on input signal level (there are no absolute thresholds) 
and does not use any type of psychoacoustics methods (e.g. masking) except onset detection.

The algorithm is calibrated for 48000 Hz, stereo. It supports any sample rate and number of channels, 
but currently sounds best on audio data with aforementioned parameters.

## Demos
**Disclamer: the following collection of files contains results based on copyrighted material, for which I own no rights.
  Short snippets of processed copyrighted material to illustrate the quality of a research algorithm on
  real-world material may be considered fair use.**

[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w). 

### Quality report (including material not in set)
- Cibo Matto — Sci-Fi Wasabi (mp3 320k): most transient detections are wrong.
- Kelela — idea 1 (mp3 320k): slap delays on vocals are smudged. This is the worst case scenario for any 
  time stretch algorithm. Elastique gives more clarity and presence.

## Notes
- `-q -1` is not necessarily should be used for speedup. 
  In many cases, faster and greedy approximation sounds better than default PGHI.
  On hard electronica, such as Carpenter Brut, it is the only option.
  Even on hiphop, such as Joji — Tick Tock, it sounds considerably better than -q 0.
- TODO: Calibrate warper's base nfft and partial derivative formulas.
- TODO: Input NaN detection and removal.
- TODO: Extract dspio/wavio to an independent package.
- A new onset detector is needed. This one is not that good.
- **Onset detection and phasor generation can be performed while warping**.
- Hypothesis: OfflineGrainReader and OfflineToOfflineGrainWriter to be removed. 
  Plain GrainReader and RegularToOfflineGrainWriter can be paired.
- We need a regression testing GitHub CI.
- Resamplers: https://codeberg.org/BillyDM/awesome-audio-dsp/src/branch/main/content/deip.pdf
- Formant shifting must be implemented after streaming.
- Momentary phase resets 
  - Are ineffective. At least if resetting at each transient ridge for the current frame.
  - Analysis lookahead will help in correct ridge detection.
  - Phase could be reset on PGHI-detected transients.
  - Phase could be reset at the start of each PGHI-detected tonal trajectory 
    (when trace\[w\] == 1).
  - We may limit amount of reset-continued ridges to, say, loudest 10-20 using 
    existing arrow data.
- Differentiation of major (full) and minor (with continued partials) phase resets.
- Discrete partial phase derivatives may perform better than reassignment.
- cmd/nanowarp: FLAC output (https://github.com/mewkiz/flac)
  - Use a clone of FabFilter Pro-L 2 “Modern” mode by [Tokyo Dawn Labs](https://t.me/tokyodawnlabsru/547) for limiting.
    Limit at 8x oversampling (requires a resampler).
    Use it by default for FLAC exports instead of clipping.
  - Take 64-bit dither from one of Airwindows plugins. Acknowledge Chris properly.
- cmd/nanowarp: allow cuts in timemap, force phase reset on each cut.
- cmd/nanowarp: Ableton Live Clip (.asd) to timemap converter.
- cmd/nanowarp: allow external onset detectors. Already possible with right timemap, 
  algorithm does phase reset on any region with Dy = 1.
- Phase ramp monotonicity is not needed. We never use `(*Curve).Sample`.

## Known issues
- BUG Writing is quantized by 1024-sample blocks. Expect loss of data at the end of file.
- No pitch modification. Requires a good resampler library,  e.g. r8brain. 
  Either port it or use through cgo.
- Slow. ≈10 seconds of output per second on Ryzen 7 7700x.
  - Even slower now because of unbufferized wavio.Encoder.GrainSeek.

## AI use disclosure
`rankfilt.go` was initially translated from C++ by free ChatGPT.
Offline modes of `GrainReader`/`GrainWriter` and some tests in `dspio` package were generated by Claude Opus (I asked a friend to do it).

Everything else was written without use of AI, if not explicitly stated otherwise.

## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)

See also https://github.com/ltfat/pvdoneright and https://github.com/y-fujii/mini_pvdr

2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Böck, S., & Widmer, G. (2013, September). Maximum filter vibrato suppression for onset detection. In Proc. of the 16th Int. Conf. on Digital Audio Effects (DAFx). Maynooth, Ireland (Sept 2013) (Vol. 7, p. 4). Citeseer.](https://www.cp.jku.at/research/papers/Boeck_Widmer_DAFx_2013.pdf)
4. [Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a real-time realization for a commercial music product.](https://thesis.unipd.it/bitstream/20.500.12608/16470/1/tesi.pdf)

## Future: Use in C/C++ through Cgo

TODO C-compatible code is not yet done. It will be in nanowarp/c subdirectory.

For future compatibility with any possible Go code, you should create your own Go library, 
and then compile it to a static C-compatible library using
```
go build -buildmode=c-archive
```

According to https://github.com/draffensperger/go-interlang/, `go build` will generate 
a header with declarations for all exported functions and a compiled object file. Use
the header and link the object, no extra system headers required.

Don't use gccgo as it is obsolete.
See also [documentation on Go interop with C](https://pkg.go.dev/cmd/cgo#hdr-C_references_to_Go).

Any Go code can not be used in an audio plugin (e.g. VST), since every new copy of a plugin (or even
any other plugin using Go code) will create it's own copy of a Go runtime, which is not supported.
Go in plugins require a runtime that supports being dynamically linked a la .NET. 
Having a bunch of go1.xx.x.dll in %WINDIR% would be very funny.

