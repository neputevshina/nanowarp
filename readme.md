# Nanowarp
An ongoing attempt to create a professional grade audio time stretching algorithm.
Reference implementation is going to be in Go, and a possible (Fil-)C implementation will share
this repo with a Go version.

Includes a modified version of github.com/youpy/go-wav (ISC license) with added 32-bit 
float WAV support export.

Current state: algorithm mostly ready. No user-facing API exists yet.

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

Nanowarp is a phase gradient heap integration (PGHI) phase vocoder[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two.

Like in original implementation of PGHI-PV, FFT is oversampled by factor of 2 with zero-padding. 
Stereo coherence is obtained through stretching mono and adding phase difference of 
respective side channels to it after stretching to get stereo signals back[4].

A phase ramp for the entire output signal is generated. Onsets are detected using rectified 
complex-domain novelty function[3]. If onset is detected, phase ramp will have a derivative of 
1 in a region around detected onset. Starting points of these sample regions are scaled by the 
stretch size, and points between regions are linearly interpolated.

Then the large-grained (nfft=4096) PGHI phase vocoder is applied, using phase ramp for the input 
sample indexes. If the derivative of the signal is 1, samples are passed through to the output 
unmodified.

## Demos
[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w)

## TODO
### Beat-emphasis onset detection

See https://www.researchgate.net/profile/Matthew-Davies-5/publication/221016733_Towards_a_musical_beat_emphasis_function/links/54465fbd0cf2d62c304db658/Towards-a-musical-beat-emphasis-function.pdf

### Testing strategy
- Various impulse train signals
- LFO FM Sine
- Vocals under hard saturation
- Drum loops
- Full tracks: pop, electronica, acoustica, black metal



## Known issues
- No pitch modification. Requires a good resampler library,  e.g. r8brain. 
  Either port it or use through cgo.
- No streaming support. All processing is in-memory with obvious RAM costs.
- Slow.
- Smaller dynamic range comparing to original audio. Related to phase accuracy.
- Phase interruption on transients sometimes results in a “chopped” sound.

## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)
2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Duxbury, C., Bello, J. P., Davies, M., & Sandler, M. (2003, September). Complex domain onset detection for musical signals. In Proc. Digital Audio Effects Workshop (DAFx) (Vol. 1, pp. 6-9). London: Queen Mary University.](https://www.dafx.de/paper-archive/2003/pdfs/dafx81.pdf)
4. [Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a real-time realization for a commercial music product.](https://thesis.unipd.it/bitstream/20.500.12608/16470/1/tesi.pdf)
