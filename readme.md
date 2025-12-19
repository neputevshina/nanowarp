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

To reduce smearing of transients, a variant of median harmonic-percussive source separation
(HPSS)[3] with a very short (nfft=512) asymmetric[5] window is first applied to the signal. Extracted impulsive components
of the signal are then warped with smaller FFT grain (64 in this case), and harmonic portion
is then warped using large grain (2048). Like in original implementation of PGHI-PV, FFT is
oversampled by factor of 2 with zero-padding. Phase is reset on princarg({d\phi^2 \over d\omega dt}) > \pi[4]
to hopefully make numerical errors from reassignment smaller. Stereo coherence is obtained through stretching 
mono and adding phase difference of respective side channels to it[6].

## Demos
[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w)

## Known issues
- No pitch modification. Requires a good resampler library,  e.g. r8brain. 
  Either port it or use through cgo.
- No streaming support. All processing is in-memory with obvious RAM costs.
- Slow. ~~Mostly from container/heap.~~ Done, now it's simply unoptimized and not vectorized.
- Preserved stereo image at the cost of metallic artifacts. 
  Must enhance existing or find another way of preserving stereo.
  Still sounds better than Elastiqué though.
- Smaller dynamic range comparing to original audio. Related to phase.
- Artifacts. Artifacts artifacts artifacts. Different after each fix.


## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)
2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.](https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)
4. [Röbel, A. (2003). A new approach to transient processing in the phase vocoder.](https://hal.science/hal-01161124/document)
5. [Olli Niemitalo's asymmetric window.](https://dsp.stackexchange.com/questions/2337/fft-with-asymmetric-windowing)
6. [Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a real-time realization for a commercial music product.](https://thesis.unipd.it/bitstream/20.500.12608/16470/1/tesi.pdf)