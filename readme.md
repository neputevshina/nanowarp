# Nanowarp
Studio-grade audio time stretching algorithm.

Reference implementation is going to be in Go, and a possible C implementation will share
this repo with a Go version.

Includes modified version of github.com/youpy/go-wav (ISC license) with added 32-bit 
float WAV support export. © 2013–2025 youpy.

Current state: algorithm done. No streaming, pitching and user-facing API exists yet.

## Installation and usage

1. Install [Go](https://go.dev/)
2. Open terminal/shell/command prompt.
3. Install Nanowarp
```
go install github.com/neputevshina/nanowarp/cmd/nanowarp@latest
```
4. Use it
```
nanowarp -p -i inputfile.wav -t <stretch> [-o outputfile.wav]
```
or
```
nanowarp -p -i inputfile.wav -from <bpm> -to <bpm> -st <semitones> [-o outputfile.wav]
```
or
```
nanowarp -p -i inputfile.wav -timemap <time map file> [-o outputfile.wav]
```
If your system can't find `nanowarp` executable, you have probably changed PATH variable in your system.
Probably the simplest way to bring it back if you are under Windows is by reinstalling the Go.
On Linux, you should know what to do.

Consult 
```
nanowarp -help
```
to get the list of available options.

## Implementation

Nanowarp is a phase gradient heap integration (PGHI) phase vocoder (aka PVDR)[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two for simplicity.

A phase ramp for the entire output signal is generated. Onsets are detected using rectified 
complex-domain novelty function[3]. If onset is detected, phase ramp will have a derivative of 
1 in a region around detected onset. Starting points of these sample regions are scaled by the 
stretch size, and points between regions are linearly interpolated.

Then the large-grained (nfft=4096) PVDR is applied, using phase ramp for input 
sample indexes. Steady portions of the signal are detected inside the same PGHI process 
by counting directions from where the phase must be integrated and their regions of influence.
If the derivative of the signal is 1, non-steady portions of the spectrum 
are bypassed to the output, and steady are integrated further.

Like in original implementation of PVDR, FFT is oversampled by factor of 2 with zero-padding. 
Stereo coherence is obtained through stretching mono and adding complex phase difference of 
respective side channels after stretching back[4].

The algorithm does not depend on input signal level (there are no absolute thresholds) 
and does not use any type of psychoacoustics methods (e.g. masking) except onset detection.

## Demos
~~[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w)~~. Obsolete.

## Notes
- Resamplers: https://codeberg.org/BillyDM/awesome-audio-dsp/src/branch/main/content/deip.pdf
- Formant shifting must be implemented after streaming.
- Phase could be reset on PGHI-detected transients.
- Phase could be reset at the start of each PGHI-detected tonal trajectory 
  (when trace\[w\] == 1).
- Analysis lookahead will help in correct ridge detection.
- We may limit amount of reset-continued ridges to, say, loudest 10-20 using 
  existing arrow data.
- Differentiation of major (full) and minor (with continued partials) phase resets.
- Discrete partial phase derivatives may perform better than reassignment.
- cmd/nanowarp: FLAC output (https://github.com/mewkiz/flac)
- cmd/nanowarp: allow cuts in timemap, force phase reset on each cut.
- cmd/nanowarp: Ableton Live Clip (.asd) to timemap converter.
- cmd/nanowarp: allow external onset detectors. Already possible with right timemap, 
  algorithm does phase reset on any region with Dy = 1.
- Cibo Matto — Sci-Fi Wasabi (mp3 320k): most transient detections are wrong.
- Phase ramp monotonicity is not needed. We never use `(*Curve).Sample`.

## Known issues
- No pitch modification. Requires a good resampler library,  e.g. r8brain. 
  Either port it or use through cgo.
- No streaming support. All processing is in-memory with obvious RAM costs.
- Slow. ≈10 seconds of output per second on Ryzen 7 7700x.
- Triple echo in time on extreme (>4x) stretches. 
  The bane of all PVDR-based algorithms due to extreme stretching of magnitude spectrum.
  Mitigated either by [SELEBI](https://arxiv.org/abs/2602.16421) or by factorization of stretch coefficient (hint from Elastiqué SDK docs).
  From `f, e := math.Frexp(stretch)`, stretch by two `e-1` times, and finish with `f*2`.
  If `e-1` is negative, shrink by `1-e` instead.
- Triple echo in frequency on high-frequency content. Can be seen on 2x stretched log sweep.
- [Modifies the tonal balance of the material.](https://mega.nz/file/emQkAArB#_HzQqUP_-1f_C9jzMcZLxSM8W21_YZoqkDXltqZgX6E) 
  Elastiqué doesn't do that.

## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)
 see also https://github.com/ltfat/pvdoneright and https://github.com/y-fujii/mini_pvdr
2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Duxbury, C., Bello, J. P., Davies, M., & Sandler, M. (2003, September). Complex domain onset detection for musical signals. In Proc. Digital Audio Effects Workshop (DAFx) (Vol. 1, pp. 6-9). London: Queen Mary University.](https://www.dafx.de/paper-archive/2003/pdfs/dafx81.pdf)
4. [Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a real-time realization for a commercial music product.](https://thesis.unipd.it/bitstream/20.500.12608/16470/1/tesi.pdf)
