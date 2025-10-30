%% Multi-Rate EEG Signal Upsampling & Alignment (CV-Ready)
% =======================================================
% - Simulates 3 EEG channels at 10, 40, 83 Hz (wearable vs. clinical)
% - Upsamples all to 83 Hz using 4 methods
% - Compares vs. 1000 Hz ground truth
% - Metrics: MSE, PSNR, spectral correlation
% - Real-world use: BCI, sleep staging, sensor fusion

clc; clear; close all;

%% 1) Ground Truth: High-Res EEG (1000 Hz)
fs_gt = 1000;
t_gt = 0:1/fs_gt:5;
N_gt = length(t_gt);

% Realistic EEG: alpha, beta, noise, artifact
alpha = 15 * sin(2*pi*10*t_gt);
beta  = 8  * sin(2*pi*20*t_gt + pi/3);
noise = 3  * randn(size(t_gt));
artifact = zeros(size(t_gt));
artifact(2000:2050) = 50;

eeg_gt = alpha + beta + noise + artifact;

%% 2) Downsample to Real Rates
fs = [10, 40, 83];
signals = cell(3,1);
timez   = cell(3,1);

for i = 1:3
    idx = round(1:fs_gt/fs(i):N_gt);
    signals{i} = eeg_gt(idx);
    timez{i}   = t_gt(idx);
end

%% 3) Plot Original
figure(1), clf, hold on;
colors = 'kmb'; markers = 'osd';
for i = 1:3
    plot(timez{i}, signals{i}, [colors(i) markers(i) '-'], ...
        'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'MarkerSize', 6);
end
xlabel('Time (s)'); ylabel('Amplitude (µV)');
title('Multi-Rate EEG (10, 40, 83 Hz)');
legend({'10 Hz', '40 Hz', '83 Hz'}, 'Location', 'northoutside');
grid on; axis tight;
axlims = axis;

%% 4) Upsample to Fastest (83 Hz)
[newSrate, fastestIdx] = max(fs);
target_N = round(length(signals{fastestIdx}) * (newSrate / fs(fastestIdx)));
newTime = (0:target_N-1) / newSrate;

%% 5) Interpolation Methods
methods = {'spline', 'linear', 'pchip', 'resample'};
nM = length(methods);
upsampled = zeros(3, target_N, nM);

for i = 1:3
    x = timez{i};
    y = signals{i};
    
    % Spline, Linear, PCHIP
    for m = 1:3
        F = griddedInterpolant(x, y, methods{m});
        upsampled(i,:,m) = F(newTime);
    end
    
    % resample()
    [p, q] = rat(newSrate / fs(i));
    sig_res = resample(y, p, q);
    t_res = (0:length(sig_res)-1)/newSrate;
    upsampled(i,:,4) = interp1(t_res, sig_res, newTime, 'linear', 'extrap');
end

%% 6) Ground Truth at 83 Hz
gt_target = interp1(t_gt, eeg_gt, newTime, 'spline');

%% 7) Spectral Analysis (FIXED: Match lengths)
figure(3), clf;

for i = 1:3
    subplot(3,1,i); hold on;
    
    % Original PSD
    [pxx_orig, f_orig] = pwelch(signals{i}, [], [], 1024, fs(i));
    plot(f_orig, 10*log10(pxx_orig), 'k-', 'LineWidth', 1.2, 'DisplayName', 'Original');
    
    % Upsampled PSD (spline)
    [pxx_up, f_up] = pwelch(upsampled(i,:,1), [], [], 1024, newSrate);
    plot(f_up, 10*log10(pxx_up), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Upsampled');
    
    % Ground Truth PSD
    [pxx_gt, f_gt] = pwelch(gt_target, [], [], 1024, newSrate);
    plot(f_gt, 10*log10(pxx_gt), 'r--', 'LineWidth', 1.5, 'DisplayName', '1000 Hz GT');
    
    title(sprintf('Channel %d: %d to %d Hz', i, fs(i), newSrate));
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    xlim([0 50]); ylim([-60 40]);
    legend('Location', 'northeast');
    grid on;
end
sgtitle('Power Spectral Density Comparison');

%% 8) Upsampled Time Domain (Spline)
figure(2), clf, hold on;
for i = 1:3
    plot(newTime, upsampled(i,:,1), [colors(i) markers(i) '-'], ...
        'LineWidth', 1.5, 'MarkerFaceColor', colors(i), 'MarkerSize', 5);
end
plot(newTime, gt_target, 'k--', 'LineWidth', 2, 'DisplayName', '1000 Hz GT');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
title('Upsampled to 83 Hz (Spline)');
legend({'10 to 83', '40 to 83', '83 Hz', '1000 Hz GT'}, 'Location', 'northoutside');
grid on; axis(axlims);

%% 9) Error Metrics
mse = zeros(3, nM);
psnr = zeros(3, nM);
corr_coef = zeros(3, nM);

for i = 1:3
    for m = 1:nM
        err = upsampled(i,:,m) - gt_target;
        mse(i,m) = mean(err.^2);
        psnr(i,m) = 10*log10(max(gt_target)^2 / mse(i,m));
        corr_coef(i,m) = corr(upsampled(i,:,m)', gt_target');
    end
end

%% 10) Print Results
fprintf('\n=== UPSAMPLING PERFORMANCE (vs 1000 Hz GT) ===\n');
fprintf('%-8s %-8s %10s %10s %10s\n', 'From', 'Method', 'MSE', 'PSNR', 'Corr');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:3
    for m = 1:nM
        fprintf('%-8s %-8s %10.2e %10.2f %10.3f\n', ...
            [num2str(fs(i)) 'Hz'], methods{m}, mse(i,m), psnr(i,m), corr_coef(i,m));
    end
end
