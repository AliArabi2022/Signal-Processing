% Advanced Signal Upsampling Demonstration
% ========================================
% This script demonstrates sophisticated upsampling techniques in MATLAB,
% including interpolation methods and the resample function. It expands on
% basic upsampling by incorporating real-world scenarios such as processing
% biomedical signals (e.g., simulated ECG data), adding noise for realism,
% comparing multiple interpolation methods, and evaluating performance
% metrics like mean squared error (MSE) and spectral analysis.
%
% Real-world applications:
% - Biomedical signal processing: Upsampling low-resolution ECG or EEG signals
%   for improved feature extraction, artifact removal, or integration with
%   high-resolution systems.
% - Audio processing: Increasing sample rates for better audio quality in
%   digital signal processing pipelines.
% - Sensor data fusion: Upsampling sparse sensor readings (e.g., from IoT
%   devices) to align with higher-frequency data streams in industrial monitoring.
%
%
% Author: Ali Arabi Bavil
% Date: October 25, 2025

clc;
clear all;
close all;

%% Section 1: Basic Example with Synthetic Data
% Original parameters
srate = 10;  % Original sampling rate (Hz)
data = [1 4 3 6 2 19];  % Sample data points
npnts = length(data);
time = (0:npnts-1)/srate;  % Time vector

% Plot original data
figure(1), clf, hold on
plot(time, data, 'ko-', 'MarkerSize', 15, 'MarkerFaceColor', 'k', 'LineWidth', 3);
title('Basic Upsampling Comparison');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Upsampling factor
upsampleFactor = 4;
newNpnts = npnts * upsampleFactor;
newTime = (0:newNpnts-1) / (upsampleFactor * srate);
newTime(newTime > time(end)) = [];  % Trim to original duration

% Actual new sampling rate
newSrateActual = 1 / mean(diff(newTime));

% Interpolation using griddedInterpolant (spline method)
F_spline = griddedInterpolant(time, data, 'spline');
updata_spline = F_spline(newTime);

plot(newTime, updata_spline, 'rs-', 'MarkerSize', 14, 'MarkerFaceColor', 'r');

% Compare with linear interpolation
F_linear = griddedInterpolant(time, data, 'linear');
updata_linear = F_linear(newTime);

plot(newTime, updata_linear, 'g*-', 'MarkerSize', 14, 'MarkerFaceColor', 'g');

% Resample function
newSrate = 40;  % Target sampling rate (4x original)
[p, q] = rat(newSrate / srate);
updata_resample = resample(data, p, q);
newTime_resample = (0:length(updata_resample)-1) / newSrate;
updata_resample(newTime_resample > time(end)) = [];
newTime_resample(newTime_resample > time(end)) = [];

plot(newTime_resample, updata_resample, 'b^-', 'MarkerSize', 14, 'MarkerFaceColor', 'b');

% Legend and limits
legend({'Original', 'Spline Interpolated', 'Linear Interpolated', 'Resampled'}, 'Location', 'Best');
xlim([0 max(time(end), newTime(end))]);

%% Section 2: Real-World Example - Simulated Noisy Biomedical Signal (ECG-like)
% Generate a simulated ECG signal with low sampling rate
% Modified to include a transient artifact to highlight method differences
fs_low = 50;  % Low sampling rate (Hz)
t_duration = 5;  % Signal duration (seconds)
t_low = 0:1/fs_low:t_duration;  % Time vector
npnts_ecg = length(t_low);

% Simulate ECG with added transient artifact
ecg_clean = 0.5 * sin(2*pi*1*t_low) + 0.3 * sin(2*pi*3*t_low + pi/2) + ...
            0.2 * sin(2*pi*5*t_low + pi);  % Base ECG
% Add a sharp transient at t=2s to test interpolation behavior
transient = zeros(size(t_low));
transient(round(length(t_low)/2.5)) = 1.5;  % Sharp spike
ecg_clean = ecg_clean + transient;

% Add realistic noise
noise_gaussian = 0.05 * randn(size(t_low));
noise_powerline = 0.1 * sin(2*pi*60*t_low);
ecg_noisy = ecg_clean + noise_gaussian + noise_powerline;

% Upsample to target rate (500 Hz)
target_fs = 500;  % High sampling rate (Hz)
upsampleFactor_ecg = target_fs / fs_low;
[p_ecg, q_ecg] = rat(upsampleFactor_ecg);
t_high = 0:1/target_fs:t_duration;
t_high(t_high > t_low(end)) = [];

% Method 1: Spline Interpolation
F_spline_ecg = griddedInterpolant(t_low, ecg_noisy, 'spline');
ecg_spline = F_spline_ecg(t_high);

% Method 2: Resample with anti-aliasing filter
ecg_resample = resample(ecg_noisy, p_ecg, q_ecg);
t_resample = 0:1/target_fs:(length(ecg_resample)-1)/target_fs;
ecg_resample(t_resample > t_low(end)) = [];
t_resample(t_resample > t_low(end)) = [];

% New Figure 2: Overlay and Zoomed-In View
figure(2), clf;

% Subplot 1: Overlay all signals
subplot(2,1,1);
plot(t_low, ecg_noisy, 'k-', 'LineWidth', 2, 'DisplayName', 'Original (Low-Sampled)');
hold on;
plot(t_high, ecg_spline, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Spline Interpolated');
plot(t_resample, ecg_resample, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Resampled');
title('Overlay of Upsampling Methods');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Location', 'Best');
grid on;

% Subplot 2: Zoomed-In View (around the transient at t=2s)
zoom_window = [1.8 2.2];  % Focus around the transient
subplot(2,1,2);
plot(t_low, ecg_noisy, 'k-', 'LineWidth', 2, 'DisplayName', 'Original (Low-Sampled)');
hold on;
plot(t_high, ecg_spline, 'r-', 'LineWidth', 1.5,'DisplayName', 'Spline Interpolated');
plot(t_resample, ecg_resample, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Resampled');
xlim(zoom_window);
title('Zoomed-In View (Near Transient Artifact)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Location', 'Best');
grid on;

%% Section 3: Performance Evaluation
% Generate ground truth: High-sampled clean ECG without noise
fs_high_gt = 1000;  % Reference sampling rate
t_gt = 0:1/fs_high_gt:t_duration;
ecg_gt = 0.5 * sin(2*pi*1*t_gt) + 0.3 * sin(2*pi*3*t_gt + pi/2) + ...
         0.2 * sin(2*pi*5*t_gt + pi);
transient_gt = zeros(size(t_gt));
transient_gt(round(length(t_gt)/2.5)) = 1.5;  % Same transient
ecg_gt = ecg_gt + transient_gt;

% Downsample ground truth to low rate
ecg_gt_down = interp1(t_gt, ecg_gt, t_low, 'linear');

% Upsample back using methods
F_spline_gt = griddedInterpolant(t_low, ecg_gt_down, 'spline');
ecg_spline_gt = F_spline_gt(t_gt(1:length(t_gt)-mod(length(t_gt), upsampleFactor_ecg)));

ecg_resample_gt = resample(ecg_gt_down, p_ecg, q_ecg);
t_resample_gt = 0:1/target_fs:(length(ecg_resample_gt)-1)/target_fs;
ecg_resample_gt_interp = interp1(t_resample_gt, ecg_resample_gt, ...
    t_gt(1:length(t_gt)-mod(length(t_gt), upsampleFactor_ecg)), 'linear');

% Compute MSE and PSNR
mse_spline = mean((ecg_gt(1:length(ecg_spline_gt)) - ecg_spline_gt).^2);
mse_resample = mean((ecg_gt(1:length(ecg_resample_gt_interp)) - ecg_resample_gt_interp).^2);

% PSNR (Peak Signal-to-Noise Ratio)
peak_signal = max(ecg_gt)^2;
psnr_spline = 10 * log10(peak_signal / mse_spline);
psnr_resample = 10 * log10(peak_signal / mse_resample);

fprintf('Performance Metrics:\n');
fprintf('MSE for Spline Interpolation: %.4f\n', mse_spline);
fprintf('MSE for Resample: %.4f\n', mse_resample);
fprintf('PSNR for Spline Interpolation: %.2f dB\n', psnr_spline);
fprintf('PSNR for Resample: %.2f dB\n', psnr_resample);

% New Figure 3: Error Analysis
figure(3), clf;
subplot(2,1,1);
error_spline = abs(ecg_spline - interp1(t_gt, ecg_gt, t_high, 'linear'));
error_resample = abs(ecg_resample - interp1(t_gt, ecg_gt, t_resample, 'linear'));
plot(t_high, error_spline, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Spline Error');
hold on;
plot(t_resample, error_resample, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Resample Error');
title('Absolute Error vs. Ground Truth');
xlabel('Time (s)');
ylabel('Absolute Error');
legend('Location', 'Best');
grid on;

% Zoomed-in error around transient
subplot(2,1,2);
plot(t_high, error_spline, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Spline Error');
hold on;
plot(t_resample, error_resample, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Resample Error');
xlim(zoom_window);
title('Zoomed-In Absolute Error (Near Transient Artifact)');
xlabel('Time (s)');
ylabel('Absolute Error');
legend('Location', 'Best');
grid on;

%% Section 4: Spectral Analysis (Enhanced)
figure(4), clf, hold on;
[pxx_orig, f_orig] = pwelch(ecg_noisy, [], [], [], fs_low);
plot(f_orig, 10*log10(pxx_orig), 'k-', 'LineWidth', 2, 'DisplayName', 'Original');

[pxx_spline, f_spline] = pwelch(ecg_spline, [], [], [], target_fs);
plot(f_spline, 10*log10(pxx_spline), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Spline');

[pxx_resample, f_resample] = pwelch(ecg_resample, [], [], [], target_fs);
plot(f_resample, 10*log10(pxx_resample), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Resample');

title('Power Spectral Density Comparison (Zoomed 0-10 Hz)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Location', 'Best');
grid on;
xlim([0 26]);  % Focus on low frequencies (ECG-relevant)
ylim([-60 0]);  % Adjust for better visibility of differences
%% Section 3: Performance Evaluation
% Compare methods using Mean Squared Error (MSE) against a high-res ground truth
% Generate ground truth: High-sampled clean ECG without noise
fs_high_gt = 1000;  % Even higher for reference
t_gt = 0:1/fs_high_gt:t_duration;
ecg_gt = 0.5 * sin(2*pi*1*t_gt) + 0.3 * sin(2*pi*3*t_gt + pi/2) + ...
         0.2 * sin(2*pi*5*t_gt + pi);

% Downsample ground truth to low rate, then upsample and compare
ecg_gt_down = interp1(t_gt, ecg_gt, t_low, 'linear');  % Simulate original sampling

% Upsample back using methods
F_spline_gt = griddedInterpolant(t_low, ecg_gt_down, 'spline');
ecg_spline_gt = F_spline_gt(t_gt(1:length(ecg_gt)-mod(length(ecg_gt), upsampleFactor_ecg)));

ecg_resample_gt = resample(ecg_gt_down, p_ecg, q_ecg);
t_resample_gt = 0:1/target_fs:(length(ecg_resample_gt)-1)/target_fs;
% Interpolate to match gt length for MSE
ecg_resample_gt_interp = interp1(t_resample_gt, ecg_resample_gt, t_gt(1:length(ecg_gt)-mod(length(ecg_gt), upsampleFactor_ecg)), 'linear');

% Compute MSE
mse_spline = mean((ecg_gt(1:length(ecg_spline_gt)) - ecg_spline_gt).^2);
mse_resample = mean((ecg_gt(1:length(ecg_resample_gt_interp)) - ecg_resample_gt_interp).^2);

fprintf('Performance Metrics:\n');
fprintf('MSE for Spline Interpolation: %.4f\n', mse_spline);
fprintf('MSE for Resample: %.4f\n', mse_resample);

%% Section 4: Spectral Analysis (Frequency Domain Comparison)
% Real-world: Useful for verifying that upsampling preserves frequency content
% without introducing artifacts, e.g., in vibration analysis for predictive maintenance.

% Compute Power Spectral Density (PSD) using Welch's method
figure(5), clf, hold on;

[pxx_orig, f_orig] = pwelch(ecg_noisy, [], [], [], fs_low);
plot(f_orig, 10*log10(pxx_orig), 'k-', 'LineWidth', 2, 'DisplayName', 'Original');

[pxx_spline, f_spline] = pwelch(ecg_spline, [], [], [], target_fs);
plot(f_spline, 10*log10(pxx_spline), 'r-', 'LineWidth', 2, 'DisplayName', 'Spline');

[pxx_resample, f_resample] = pwelch(ecg_resample, [], [], [], target_fs);
plot(f_resample, 10*log10(pxx_resample), 'b-', 'LineWidth', 2, 'DisplayName', 'Resample');

title('Power Spectral Density Comparison');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Location', 'Best');
grid on;
xlim([0 150]);  % Focus on relevant frequencies for ECG

% End of script