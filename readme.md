# Nanowarp
An ongoing attempt of professional grade audio time stretching algorithm.
Reference implementation is going to be in Go, and a possible C implementation will share
this repo with a Go version.

Includes a modified version of github.com/youpy/go-wav (ISC license) with added 32-bit 
float WAV support export.

## Implementation
Nanowarp is a phase gradient heap integration (PGHI) phase vocoder[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two.

To reduce smearing of transients, a variant of median harmonic-percussive source separation
(HPSS)[3] with a very short window is first applied to the signal. Extracted impulsive components
of the signal are then warped with smaller FFT grain (512 in this case), and harmonic portion
is then warped using large grain (4096). Like in original implementation of PGHI-PV, FFT is
oversampled by factor of 2 with zero-padding.

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
3. [Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.](https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)

