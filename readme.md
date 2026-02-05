# Nanowarp
An ongoing attempt to create a professional grade audio time stretching algorithm.
Reference implementation is going to be in Go, and a possible (Fil-)C implementation will share
this repo with a Go version.

Includes a modified version of github.com/youpy/go-wav (ISC license) with added 32-bit 
float WAV support export.

Current state: not ready. Algorithm is still in the state of polishing, user-facing API does not exist.

## Installation and usage

```
$ go install github.com/neputevshina/nanowarp/cmd/nanowarp@latest
$ nanowarp -i inputfile.wav -t <stretch> [-o outputfile.wav]
```

## Implementation

Nanowarp is a phase gradient heap integration (PGHI) phase vocoder[1] where partial derivatives 
of phase are obtained through time-frequency reassignment[2]. This way accurate phase-time 
advance can be obtained using only one windowed grain instead of two.

Like in original implementation of PGHI-PV, FFT is oversampled by factor of 2 with zero-padding. 
Stereo coherence is obtained through stretching mono and adding phase difference of 
respective side channels to it after stretching to get stereo signals back[6].

### -alt=false
To reduce smearing of transients, a variant of median harmonic-percussive source separation
(HPSS)[3] with a very short (nfft=512) asymmetric[5] window is first applied to the signal. 
Extracted impulsive components of the signal are then warped with smaller FFT grain (64 in 
this case), and harmonic portion is then warped using large grain (4096). 
Phase is fully reset on onset points, detected from the percussive signal separately.

### -alt=true
Median HPSS with same asymmetric window and is used to extract percussive component of the signal.
HPSS mask is eroded by frequency, and some certain conditions applied to exclude most of the 
incorrectly classified noise.

Then, a phase ramp for the entire output signal is generated. If percussive signal is non-zero in 
the time domain, phase ramp will have a derivative of 1 in these samples. A switch can be from 15 to 
50 ms long max. Indexes of percussive points are scaled by the stretch size, points between switches
are linearly interpolated.

Then the same large PGHI phase vocoder is applied, using phase ramp for the input sample indexes.
If the derivative of the signal is 1, samples are passed through to the output.

> TODO: Lite mode: reset the phase at peaks of cross-correlation of analyzed and synthesized signals.

## Demos
~~[Listen here](https://mega.nz/folder/ayZwxaAA#pcw2-oE-lwXRmPC6g4fg6w)~~ Outdated

## TODO Testing strategy
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
- Amplidude modulation on transients. A consequence of the previous issue.
- Artifacts. Artifacts artifacts artifacts. Different after each fix.
- `-alt=true`: phase interruption on transients results in a “chopped” sound.


## References
1. [Průša, Z., & Holighaus, N. (2017). Phase vocoder done right.](https://ltfat.org/notes/ltfatnote050.pdf)
2. [Flandrin, P. et al. (2002). Time-frequency reassignment: from principles to algorithms.](https://hal.science/hal-00414583/document)
3. [Fitzgerald, D. (2010). Harmonic/percussive separation using median filtering.](https://dafx10.iem.at/proceedings/papers/DerryFitzGerald_DAFx10_P15.pdf)
4. [Röbel, A. (2003). A new approach to transient processing in the phase vocoder.](https://hal.science/hal-01161124/document)
5. [Olli Niemitalo's asymmetric window.](https://dsp.stackexchange.com/questions/2337/fft-with-asymmetric-windowing)
6. [Altoè, A. (2012). A transient-preserving audio time-stretching algorithm and a real-time realization for a commercial music product.](https://thesis.unipd.it/bitstream/20.500.12608/16470/1/tesi.pdf)
