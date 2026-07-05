# Nanowarp
Studio-grade audio time stretching algorithm.

Reference implementation is going to be in Go, and a possible C implementation will share
this repo with a Go version.

Includes a modified version of github.com/youpy/go-wav (ISC license) with added 32-bit 
float WAV support export. © 2013–2025 youpy.

Current state: perfecting the algorithm. No user-facing API exists yet.

## Installation and usage

1. Install [Go](https://go.dev/)
2. Open terminal/shell/command prompt.
3. Install Nanowarp
```
go install github.com/neputevshina/nanowarp/cmd/nanowarp@latest
```
4. Use it
```
nanowarp -i inputfile.wav -t <stretch> [-o outputfile.wav]
```
or
```
nanowarp -i inputfile.wav -from <bpm> -to <bpm> -st <semitones> [-o outputfile.wav]
```
If your system can't find `nanowarp` executable, you have probably changed PATH variable in your system.
Probably the simplest way to bring it back if you are under Windows is by reinstalling the Go.
On Linux, you should probably know what to do.

Consult 
```
nanowarp -help
```
to get the list of available options.

## Implementation

Nanowarp is a phase gradient heap integration (PGHI) phase vocoder (aka PVDR)[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two for simplicity.

Like in original implementation of PVDR, FFT is oversampled by factor of 2 with zero-padding. 
Stereo coherence is obtained through stretching mono and adding complex phase difference of 
respective side channels after stretching back[4].

A phase ramp for the entire output signal is generated. Onsets are detected using rectified 
complex-domain novelty function[3]. If onset is detected, phase ramp will have a derivative of 
1 in a region around detected onset. Starting points of these sample regions are scaled by the 
stretch size, and points between regions are linearly interpolated.

Then the large-grained (nfft=4096) PVDR is applied, using phase ramp for the input 
sample indexes. If the derivative of the signal is 1, samples are passed through to the output 
*essentially* unmodified.

The algorithm does not depend on input signal level (there are no absolute thresholds) 
and does not use any type of psychoacoustics methods (e.g. masking) except onset detection.

## Demos
~~[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w)~~. Obsolete.

## Notes
- There could be ways to optimize the “phase gradient heap integration” to not need a heap.
- Some more onset detectors:
  - https://www.cp.jku.at/research/papers/Boeck_Widmer_DAFx_2013.pdf
  - https://www.dlsi.ua.es/~pertusa/pub/pdf/ciarp05.pdf
  - Expecting a regular beat might be bad for some types of music.
- SELEBI exists (preprint): https://arxiv.org/abs/2602.16421
- ~~PGHI, being a “brute-force sinusoidal modeling”, probably can be abused as a tonality measure for ruling out erroneous onset detections.~~ 
  ~~It can't, but it's still a cool concept to keep in mind.~~
  It actually can, but you need to see where points are coming from, not where they lead.
- [Non-causal PGHI](https://ltfat.org/notes/ltfatnote040.pdf) is ineffective because PGHI integrates the phase locally, 
  ignoring overlap, so it is impossible to obtain globally coherent phase with phase resets using this method. 
  We need some way to use the phase of up to overlap number of frames.
- From the cellular automaton/WFC view, increasing the neighborhood of PGHI is ineffective. 
  Anything other than current von Neumann (including Moore) neighborhood gives worse results both in causal and non-causal modes.
  Can be useful for more elegant theoretical definition of PGHI.
- Short-time (3×olap frames) phase reconstruction (Griffin-Lim and friends) is ineffective for eliminating clicks.
  I don't know how, but it only makes worse. Skill issue maybe.
- Resamplers: https://codeberg.org/BillyDM/awesome-audio-dsp/src/branch/main/content/deip.pdf
- Formant shifting must be implemented after streaming.

## Known issues
- No pitch modification. Requires a good resampler library,  e.g. r8brain. 
  Either port it or use through cgo.
- No streaming support. All processing is in-memory with obvious RAM costs.
- Slow. ≈10 seconds of output per second on Ryzen 7 7700x.
- Interruptions because of phase reset at incorrect onset detections.
  May be fixed by simply not resetting at onsets, but it significantly impairs output quality.
- Does not reconstruct the signal perfectly,
  DC turns into a slow oscillation after FFT/IFFT cycle and is not equal
  to doubly applied windowing.
  Wrong Blackman-Harris window usage (doesn't occur with Hann IIRC)
  or a bug in gonum/fourier (unlikely).
- Triple echo in time on extreme (>4x) stretches. 
  The bane of all PVDR-based algorithms because of extreme stretching of the magnitude spectrum.
  Mitigated by factorization of stretch coefficient and repeated stretching (hint from Elastiqué SDK docs).
- Triple echo in frequency on high-frequency content. Can be seen on 2x stretched log sweep.
- [Modifies the tonal balance of the material.](https://mega.nz/file/emQkAArB#_HzQqUP_-1f_C9jzMcZLxSM8W21_YZoqkDXltqZgX6E) 
  Elastiqué doesn't do that.

## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)
 see also https://github.com/ltfat/pvdoneright and https://github.com/y-fujii/mini_pvdr
2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Duxbury, C., Bello, J. P., Davies, M., & Sandler, M. (2003, September). Complex domain onset detection for musical signals. In Proc. Digital Audio Effects Workshop (DAFx) (Vol. 1, pp. 6-9). London: Queen Mary University.](https://www.dafx.de/paper-archive/2003/pdfs/dafx81.pdf)
4. [Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a real-time realization for a commercial music product.](https://thesis.unipd.it/bitstream/20.500.12608/16470/1/tesi.pdf)
