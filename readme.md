# Nanowarp
An ongoing attempt of professional grade audio time stretching algorithm.
Reference implementation is going to be in Go, and a possible C implementation will share
this repo with a Go version.

## Implementation
Nanowarp is a phase gradient heap integration (PGHI) phase vocoder[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two. Horizontal phase 
advance in this version of PGHI is allowed only for bins with local maxima in magnitude,
working as a simplest form of auditory noise masking method.

To reduce smearing of transients, a variant of median harmonic-percussive source separation
(HPSS) with a very short window is first applied to the signal. Extracted impulsive components
of the signal are then warped with smaller FFT grain (512 in this case), and harmonic portion
is then warped using large grain (4096). Like in original implementation of PGHI-PV, FFT is
oversampled by factor of 2 with zero-padding.

**Currently none of filter sizes (FFT and HPSS quantile filters) are corrected for sample rate.**
**Current implementation is calibrated for 48kHz.**

## Demos
- [Joji — Tick Tock (original)](https://mega.nz/file/m7AGCbCR#FKddyqQFEPtKz6sTvZcs1mz4FBbrrgOlAePzoxdHWpE)
- [Joji — Tick Tock (2x stretch)](https://mega.nz/file/v3RADSLD#wkVLz_feIvGGKMrE_jAeKxhsjbyRsSIXW09wKiXn3_A)
- [Joji — Tick Tock (0.5x stretch)](https://mega.nz/file/Kvx03RAI#4wxdcfG0xHNtPPcgzDjJv2VGoDYYrC4wx3wZYTKw89Q)
- ~~Joji — Tick Tock (+12 pitch, external)~~
- [Boa — Welcome (original)](https://mega.nz/file/26BSAAyL#MMiGQ_cGTgW568C9CRbGi0mSpK3eLoOXLqQ03SYhyDQ)
- [Boa — Welcome (2x stretch)](https://mega.nz/file/rnhGTR6Y#AiwE8upG0ETx1rTWbKUoBrY0Hi8g1viUpucibYYZABY)
- [Boa — Welcome (0.5x stretch)](https://mega.nz/file/f7522Qba#UKvFU1lFHDyloxI-5j2rYMbT1HHgQovVC2XVlY7aC48)
- ~~Boa — Welcome (-12 pitch, external)~~

## References
[1] Průša, Z., & Holighaus, N. (2017, August). Phase vocoder done right. In 2017 25th European
    signal processing conference (Eusipco) (pp. 976-980). IEEE.
[2] Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.
[3] Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.
