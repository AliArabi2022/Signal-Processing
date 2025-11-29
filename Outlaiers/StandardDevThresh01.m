%% Robust Outlier Detection and Correction in Log-Normal Noisy Signals
% This script demonstrates a practical and widely used approach for cleaning
% time-series data that is mostly log-normally distributed but contaminated
% by sporadic large outliers (e.g. sensor glitches, transmission errors,
% financial price spikes, etc.).
%
% Method:
% 1. Generate a synthetic positive log-normal 
%   signal (simulating real-world multiplicative noise).
% 2. Inject a controlled number of artificial extreme
%   outliers, including both very large positive spikes
%   and negative (or near-zero) spikes, to mimic real 
%   sensor glitches, data corruption, or transmission 
%   errors.
% 3. Automatically detect outliers in both directions
%   using robust statistics:
%        - Compute the median and Median Absolute Deviation
%    (MAD) of the signal
%        - Define symmetric thresholds as:
%            lower threshold = median − k × MAD
%            upper threshold = median + k × MAD
%        (with k ≈ 5–7; far more reliable than mean ± σ 
%   on skewed or heavy-tailed data)
% 4. Flag all points lying outside these lower and upper
%   thresholds as outliers.
% 5. Replace only the flagged outlier values with linearly
%   interpolated values derived from the surrounding valid 
%   (non-outlier) points using MATLAB’s griddedInterpolant
%   (preserving local trends and continuity).
% 6. Visualize the original contaminated signal together
%   with both thresholds and the final cleaned signal for 
%   immediate quality assessment.
% 
% Author: Ali Arabi Bavil
% Date:   November 2025

clear; close all; clc;

%% Parameters (same as your original spirit)
N           = 2000;
nOutliers   = 50;                    % Number of artificial outliers
outlierMult = 10;                    % How extreme the outliers are

%% 1. Generate log-normal signal
time   = (1:N)' / N;                         % Normalized time [0–1]
signal = exp(0.5 * randn(N,1));               % Pure log-normal noise (always positive)

%% 2. Inject random outliers — both high positive and negative spikes
rng(123);  % For reproducibility
randpnts = randi(N, [nOutliers, 1]);

% Some outliers very high, some negative (like sensor going to zero or below)
outlierValues = [rand(floor(nOutliers/2),1) * range(signal) * outlierMult + max(signal);
                 -rand(ceil(nOutliers/2),1) * range(signal) * outlierMult];

signal(randpnts) = outlierValues(randperm(nOutliers));  % Shuffle them

%% 3. Plot original contaminated signal
figure(1); clf; hold on;
plot(time, signal, 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'LineWidth', 0.8);
xlabel('Normalized Time'); ylabel('Amplitude');
title('Log-Normal Signal with Positive and Negative Outliers');

%% 4. Robust bidirectional outlier detection using Median Absolute Deviation (MAD)
% This is the gold-standard simple method for skewed/heavy-tailed data

medianVal = median(signal);
madVal    = mad(signal, 1);           % MAD with normal consistency (1.4826 factor built-in)
k         = 5;                        % Common choice: 5–7 for very aggressive cleaning

lowerThreshold = medianVal - k * madVal;
upperThreshold = medianVal + k * madVal;

% Plot both thresholds
plot(time([1 end]), [1 1]*upperThreshold, 'r--', 'LineWidth', 1.5);
plot(time([1 end]), [1 1]*lowerThreshold, 'r--', 'LineWidth', 1.5);

text(0.02, upperThreshold, ' Upper threshold ', 'Color','r', 'VerticalAlignment','bottom');
text(0.02, lowerThreshold, ' Lower threshold ', 'Color','r', 'VerticalAlignment','top');

%% 5. Detect outliers on BOTH sides
outliers = signal > upperThreshold | signal < lowerThreshold;

fprintf('%d outliers detected (%.2f%% of data)\n', sum(outliers), 100*mean(outliers));

%% 6. Interpolate only the bad points — everything else unchanged!
% This block is almost identical to your original code
F = griddedInterpolant(time(~outliers), signal(~outliers), 'linear', 'nearest');

signalClean = signal;                              % Keep original
signalClean(outliers) = F(time(outliers));         % Replace only outliers

%% 7. Plot cleaned result (same style as yours)
plot(time, signalClean, 'ro-', 'LineWidth', 1.2, 'MarkerFaceColor', 'w');

legend({'Original + outliers', 'Thresholds', 'Cleaned signal'}, 'Location', 'best');
grid on;
title(sprintf('Outlier Correction via Interpolation\n(%d outliers replaced)', sum(outliers)));
