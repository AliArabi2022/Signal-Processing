clc; clear all; close all;
% Clear command window, clear workspace, close all figures.
% Good practice at the start of a script.

% number of points in signal
n = 10000;
% Total length of the synthetic signal = 10 000 samples.

% create signal
[origsig,signal] = deal( cumsum( randn(n,1) ) );
% randn(n,1)       → Gaussian white noise (zero mean, variance 1)
% cumsum(...)      → cumulative sum → creates a 1/f-like random walk
%                     (low-frequency drift + higher-frequency fluctuations)
% deal(A,A)        → copies the same vector into two variables
% origsig          → untouched original signal (kept for plotting)
% signal           → working copy that will have a gap introduced

% remove a specified window
boundaryPnts = [ 5000 6500 ];
% Define the indices where the gap starts and ends:
% Gap runs from sample 5000 to 6500 inclusive (1501 points removed)

signal(boundaryPnts(1):boundaryPnts(2)) = 0/0;
% 0/0 produces NaN in MATLAB.
% This turns the entire gap into NaN values (a common way to mark missing data)

% FFTs of pre- and post-window data
fftPre = fft(signal( boundaryPnts(1)-diff(boundaryPnts)-1 : boundaryPnts(1)-1 ));
% diff(boundaryPnts) = 1500 → length of the gap - 1
% So we take 1501 points BEFORE the gap:
%   from index 5000 - 1500 - 1 = 3499 to index 4999
% These are the last 1501 valid points before the gap.
% fftPre = Fourier transform of this pre-gap segment.

fftPst = fft(signal( boundaryPnts(2)+1 : boundaryPnts(2)+diff(boundaryPnts)+1 ));
% Take 1501 points AFTER the gap:
%   from index 6501 to 6501 + 1500 = 8001
% fftPst = Fourier transform of this post-gap segment.

% interpolated signal is a combination of mixed FFTs and straight line
mixeddata = detrend(ifft( ( fftPre + fftPst ) / 2 ));
% (fftPre + fftPst)/2 → average the two spectra (same length, same frequencies)
% ifft(...)          → inverse FFT → gives a periodic signal with average spectrum
% detrend(...)       → removes any linear trend from this periodic reconstruction
% mixeddata now contains a "spectral interpolation" of length 1501

linedata = linspace(0,1,diff(boundaryPnts)+1)' * ...
           (signal(boundaryPnts(2)+1) - signal(boundaryPnts(1)-1)) + ...
           signal(boundaryPnts(1)-1);
% Create a straight line that connects the last valid point before the gap
% (signal(4999)) to the first valid point after the gap (signal(6501)).
% linspace(0,1,N)' → column vector from 0 to 1 with N points (1501 here)
% Multiplies by the total jump in amplitude → linear ramp
% Then adds the starting value → perfect straight-line connection.
% This is classic linear interpolation across the gap.

% sum together for final result
linterp = mixeddata + linedata;
% Combine the two approaches:
%   - linedata   ensures the interpolation starts and ends at the correct values
%                (preserves level and avoids jumps)
%   - mixeddata  adds back the average spectral content (oscillatory patterns,
%                including ~7 Hz alpha-like fluctuations)
% Result: smooth transition + preserved frequency content → very good for EEG!

% put the interpolated piece into the signal
filtsig = signal;                                    % copy the gappy signal
filtsig(boundaryPnts(1):boundaryPnts(2)) = linterp;  % replace gap with interpolation
% filtsig is now the final repaired signal (NaNs are gone)

% -------------------------- Plotting --------------------------
figure(1), clf
plot(1:n, [origsig signal+5 filtsig+10])
% Plot three versions vertically offset for easy comparison:
%   - original clean random walk
%   - signal with big NaN gap (appears as blank)
%   - repaired signal (filtsig)

legend({'Original'; 'With gap'; 'Filtered'})
title('EEG-like signal with gap removal using spectral + linear interpolation')
zoom on