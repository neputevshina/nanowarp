# Nanowarp
An ongoing attempt of professional grade audio time stretching algorithm.
Reference implementation is going to be in Go, and a possible C implementation will share
this repo with a Go version.

Includes a modified version of github.com/youpy/go-wav (ISC license) with added 32-bit 
float WAV support export.

## Installation and usage

```
$ go install github.com/neputevshina/nanowarp/cmd/nanowarp
$ nanowarp -i inputfile.wav -t <stretch> [-o outputfile.wav]
```

## Implementation
Nanowarp is a phase gradient heap integration (PGHI) phase vocoder[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two.

To reduce smearing of transients, phase is resetted for impulsive parts of the spectrum[3] and 
phase reset points are stretched with smaller window size (64 samples in time). 

## Demos
[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w)

## Known issues
- No pitch modification. Requires a good resampler library,  e.g. r8brain. 
  Either port it or use through cgo.
- No streaming support. All processing is in-memory with obvious RAM costs.
- Slow. Mostly from container/heap.
- Bubbling artifacts and smeared transients. 

## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)
2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Röbel, A. (2003). A new approach to transient processing in the phase vocoder.](https://hal.science/hal-01161124/document)
4. [Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.](https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)

