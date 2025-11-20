# Signal-Processing

A structured collection of signal processing algorithms, simulations, and tutorials implemented in **MATLAB** and **Python**. The repository covers spectral analysis, wavelets, filtering, resampling, time-series denoising, and related foundational techniques used in modern communication systems and electronic signal workflows.

---

## Key Topics

- **Spectral Analysis**  
  FFT, windowing, PSD estimation, spectral leakage analysis.

- **Wavelet Processing**  
  Discrete wavelet transforms, denoising, multi-resolution decomposition.

- **Time-Series Denoising**  
  Reduction of noise using wavelets, filtering, and transform methods.

- **Convolution & Filtering**  
  FIR/IIR filter implementation, kernel-based convolution, smoothing.

- **Resampling & Interpolation**  
  Downsampling, upsampling, anti-aliasing, interpolation methods.

- **Complex Signal Processing**  
  I/Q signals, analytic representations, Hilbert transforms.

---

## Repository Structure

Signal-Processing/<br />
├── MATLAB/<br />
│ ├── Spectral/<br />
│ ├── Wavelet/<br />
│ ├── TimeSeriesDenoising/<br />
│ ├── Convolution/<br />
│ └── Resample/<br />
│<br />
├── Python/<br />
│ ├── spectral.ipynb<br />
│ ├── wavelet.ipynb<br />
│ ├── convolution.ipynb<br />
│ ├── denoising.ipynb<br />
│ └── resampling.ipynb<br />
│
├── data/<br />
├── docs/<br />
└── README.md<br />

Each subdirectory contains modular examples illustrating the theory and implementation of each technique.

---

## Getting Started

### Requirements

**MATLAB:**
- Any recent version (R2020+ recommended)

**Python:**
Install via:
```bash
pip install numpy scipy matplotlib pywt jupyter
```
Installation
```bash
git clone https://github.com/AliArabi2022/Signal-Processing.git
cd Signal-Processing
```

MATLAB: Open .m files directly.

Python: Run notebooks via:
```bash
jupyter notebook
```

## Usage

Example workflows:

Frequency Domain Analysis
Apply window functions → compute FFT → evaluate spectrum.

Wavelet Denoising
Perform DWT → thresholding → reconstruct clean signal.

Resampling
Downsample/upsample → compare aliasing → apply anti-aliasing filters.

Filtering
Design FIR/IIR filters → apply convolution → analyze stability.

Each folder contains comments, equations, and plots to make the concepts accessible.
--- 

## Roadmap

Planned additions:

Filter design toolbox (MATLAB + Python)

Adaptive filtering (LMS, NLMS, RLS)

Real-time signal streaming examples

Kalman filtering and state-space signal processing

Machine learning + signal processing examples

Expanded documentation with derivations and problem sets

## Contributing

See CONTRIBUTING.md for guidelines on pull requests, style, and code organization.

## License

This project is released under the MIT License. See the LICENSE file for details.

## Maintainer

Ali Arabi
Researcher in Communication Systems & Signal Processing
GitHub: AliArabi2022

---

