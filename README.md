# FM Signal Processing Pipeline using RTL-SDR

A complete Digital Signal Processing (DSP) pipeline for FM radio reception using Software Defined Radio (RTL-SDR). The project captures raw IQ samples, performs FM demodulation, compares FIR and IIR filtering techniques, applies adaptive Wiener filtering for noise reduction, and visualizes the signal using Short-Time Fourier Transform (STFT).

---

## Features

- Wideband FM spectrum scanning
- Automatic station selection
- IQ data acquisition using RTL-SDR
- FM demodulation
- FIR Low-Pass Filtering (Hamming Window)
- IIR Low-Pass Filtering (Chebyshev Type II)
- Wiener adaptive noise reduction
- STFT-based time-frequency analysis
- Audio reconstruction
- Performance comparison using SNR, phase response, and group delay

---

## Processing Pipeline

```
RTL-SDR
   │
   ▼
Wideband FM Scan
   │
   ▼
Station Selection
   │
   ▼
IQ Data Capture
   │
   ▼
FM Demodulation
   │
   ▼
Baseband Filtering
 ┌──────────────┐
 │              │
 ▼              ▼
FIR Filter    IIR Filter
 │              │
 └──────┬───────┘
        ▼
 Wiener Denoising
        │
        ▼
 STFT Analysis
        │
        ▼
 Audio Output
```

---

## Project Tasks

### 1. Wideband FM Scanning

- Scan the FM broadcast spectrum
- Detect the strongest broadcast stations
- Automatically select a station based on the assigned index

---

### 2. IQ Data Acquisition

- Tune RTL-SDR to the selected FM station
- Capture raw IQ samples
- Perform baseline FM demodulation

---

### 3. FIR vs IIR Filtering

Implemented two different low-pass filters for baseband processing.

#### FIR Filter

- Hamming Window
- Order 101
- Linear phase
- Constant group delay

#### IIR Filter

- Chebyshev Type II
- Order 6
- Stopband attenuation: 60 dB
- Computationally efficient
- Non-linear phase response

---

### 4. Wiener Noise Reduction

Implemented adaptive Wiener filtering using estimated signal and noise power spectra.

Advantages:

- Reduces background noise
- Improves audio clarity
- Adaptive filtering based on signal statistics

---

### 5. Time-Frequency Analysis

Applied Short-Time Fourier Transform (STFT) to visualize

- Raw FM audio
- FIR filtered signal
- Wiener filtered signal

---

## Technologies Used

- MATLAB
- RTL-SDR
- DSP Toolbox
- FFT
- STFT
- FIR Filter Design
- IIR Filter Design
- Wiener Filtering

---

## Results

The project demonstrates

- Successful FM signal acquisition
- Accurate FM demodulation
- Better signal fidelity using FIR filtering
- Reduced computational complexity using IIR filtering
- Improved audio quality through adaptive Wiener denoising
- Time-frequency visualization using spectrograms

---

## Key Concepts

- Digital Signal Processing
- Software Defined Radio
- FM Demodulation
- FIR Filters
- IIR Filters
- Wiener Filtering
- Spectrogram Analysis
- Frequency Domain Processing
- Signal-to-Noise Ratio (SNR)

---

## Future Work

- Stereo FM decoding
- RDS (Radio Data System) decoding
- Real-time SDR processing
- Adaptive filter optimization
- Automatic station classification

---

## Acknowledgements

Developed as part of the Digital Signal Processing coursework at IIIT Bangalore.
