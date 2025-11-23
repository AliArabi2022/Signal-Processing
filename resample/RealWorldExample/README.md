# DTW Speech Alignment – Studio Quality Voice Time-Warping

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-blue.svg)](https://www.mathworks.com)
[![Audio](https://img.shields.io/badge/Audio-16kHz%20Speech-green)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

https://github.com/AliArabi2022/Signal-Processing/tree/master/resample/RealWorldExample/hello_warped_STUDIO_PERFECT.MP4

**Say "helloooooo" slowly → then "hello" fast → hear them become perfectly synchronized.**  
This project demonstrates **Dynamic Time Warping (DTW)** applied to real speech with **professional-grade audio quality** — no beeps, no robot voice, natural timbre preserved.

Used in:
- Voice conversion
- Singing alignment
- Speech rate normalization
- Audio post-production
- Research demos & education

## Demo (What You’ll Hear)

| Original Slow | → | Warped Fast (aligned) |
|---------------|---|-------------------------|
| "helllllloooooooo" | → | "helllllloooooooo" (but spoken fast originally!) |

They sound **almost identical** — that’s the power of high-quality DTW warping.

## Features

- Live microphone recording (slow + fast "hello")
- Short-time energy envelope extraction
- Classic **Sakoe-Chiba DTW** with manual backtracking
- **Studio-quality time warping** using:
  - PCHIP interpolation
  - Anti-aliased filtering
  - LPC-based **formant preservation**
  - RMS loudness matching
  - Fade in/out
- 3D cost surface visualization with green optimal path
- Compares manual DTW vs MATLAB’s `dtw()`
- Saves warped audio as 24-bit WAV

## Requirements

- MATLAB R2020b or later
- Signal Processing Toolbox (for `lpc`, `butter`, `filtfilt`)
- Audio Toolbox (for `audiorecorder`, `audiowrite`)

> Works perfectly on Windows, macOS, and Linux.

## How to Run

1. Clone the repo:
```bash
git clone https://github.com/AliArabi2022/Signal-Processing/tree/master/resample/RealWorldExample/
cd RealWorldExample
```
2. Open main_DTW_Speech_Alignment.m in MATLAB
3. Run the script
4. When prompted:
    - Say "helloooooo" slowly (3 seconds)
    - Then say "hello" quickly

5. Listen to the magic!

Output: hello_warped_STUDIO.wav 

## File Structure

DTW-Speech-Alignment-MATLAB/<br>
├── main_DTW_Speech_Alignment.m     Main script (run this!)
├── warpSpeechDTW.m                 Studio-quality warping<br> function<br>
├── hello_warped_STUDIO.wav  Example output<br>
├── README.md<br>
└── LICENSE



## Example Output
After running:<br>
matlab<br>
=== DTW Summary ===<br>
Manual total cost      : 8.247813<br>
MATLAB dtw distance    : 8.247813<br>
Warping path length    : 487 points<br>
Slow utterance length  : 2.156 seconds<br>
Fast utterance length  : 0.912 seconds<br>
Speed ratio (fast/slow): 2.36×<br>