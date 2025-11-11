%% =====================================================================
%  REAL-WORLD EXAMPLE: Irregularly sampled EEG (7 Hz alpha rhythm)
%  ---------------------------------------------------------------
%  1. Simulate a clean 2-second EEG trace @ 1324 Hz containing a
%     Gaussian-shaped 7 Hz alpha peak.
%  2. Randomly drop samples according to an exponential inter-sample
%     distribution (mimics bursty packet loss).
%  3. Reconstruct a regular 1324 Hz signal with spline interpolation.
%  4. Compare time-domain traces and power spectra.
% ======================================================================
clear; close all; clc;

%% ----------------------- 0. Load .mat --------------------
File = load('subject_01.mat');          % <-- put file in MATLAB path
targetSrate = 1324;                     % Hz

% ----- USER: MANUAL PEAK FREQUENCY (change any time) -----
peakfreq      = 8;                     % <--- SET YOUR DESIRED FREQUENCY
useManualPeak = false;                   % true = manual, false = auto-detect
alphaIdx = 9;
% ---------------------------------------------------------

duration = 120;                         % seconds (first 2 min)

%% ----------------------- 1. Extract Oz (channel 16) --------------------
timeRaw = File.SIGNAL(:,1) / 1e6;       % µs → s
eegRaw  = detrend(File.SIGNAL(:,16));   % Oz channel, remove DC
% ----------------------- low pass filter --------------------
fcutoff = 14;
transw  = .2;
order   = round( 17*512/fcutoff );
shape   = [ 1 1 0 0 ];
frex    = [ 0 fcutoff fcutoff+fcutoff*transw 512/2 ] / (512/2);
% filter kernel
filtkern = firls(order,frex,shape);
eegRaw = filtfilt(filtkern,1,eegRaw);







% Trim to first `duration` seconds
keepIdx = timeRaw <= duration;
timeRaw = timeRaw(keepIdx);
eegRaw  = eegRaw(keepIdx);

fprintf('Loaded %.1f s of real EEG (Oz) @ 512 Hz\n', duration);

%% ----------------------- 2. Upsample to 1324 Hz --------------------
% Regular grid that **exactly** spans the raw interval
timevec = linspace(timeRaw(1), timeRaw(end), round(targetSrate*duration));

F_up = griddedInterpolant(timeRaw, eegRaw, 'pchip');
signal = F_up(timevec);
signal = signal(:);

%% ----------------------- 3. Irregular Sampling --------------------
rng(42);
interSample = ceil(exp(4*rand(numel(signal),1)));
sampIdx = cumsum([1; interSample]);
sampIdx(sampIdx > numel(signal)) = [];

irregSignal = signal(sampIdx);
irregTime   = timevec(sampIdx);

[irregTime, sortIdx] = sort(irregTime);
irregSignal = irregSignal(sortIdx);

%% ----------------------- 4. Spline Reconstruction --------------------
F = griddedInterpolant(irregTime, irregSignal, 'spline');
newsignal = F(timevec);
newsignal = newsignal(:);

%% ----------------------- 5. Peak Frequency (Manual / Auto) --------------------
winLen = targetSrate;
[pxx, f] = pwelch(signal, hanning(winLen), winLen/2, [], targetSrate);

if ~useManualPeak
    % ---- AUTO DETECTION (8-12 Hz) ----
    alphaIdx = f >= 8 & f <= 12;
    [~, idx] = max(pxx(alphaIdx));
    peakfreq = f(alphaIdx); peakfreq = peakfreq(idx);
    fprintf('\n>>> AUTO-DETECTED ALPHA PEAK: %.2f Hz <<<\n', peakfreq);
else
    fprintf('\n>>> USING MANUAL PEAK FREQUENCY: %.2f Hz <<<\n', peakfreq);
end

% ---- FIND CLOSEST FREQUENCY BIN (robust) ----
[~, binIdx] = min(abs(f - peakfreq));
peakfreqBin = f(binIdx);          % exact frequency that exists in PSD
peakPower   = pxx(binIdx);        % power at that bin

%% ----------------------- 6. Time-Domain Plot --------------------
figure('Color','w','Position',[100 100 1000 600]); clf; hold on;
plot(timevec, signal,      'k', 'LineWidth',1.2);
plot(irregTime, irregSignal, 'ro', 'MarkerFaceColor','r','MarkerSize',5);
plot(timevec, newsignal,   'm.', 'MarkerSize',3);
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title(sprintf('Real EEG – Peak @ %.2f Hz (%s)', peakfreq, ...
    iif(useManualPeak,'Manual','Auto')));
legend({'Real (1324 Hz)','Irregular','Reconstructed'},'Location','southoutside');
grid on; axis tight;

%% ----------------------- 7. PSD with Marker (FIXED) --------------------
avgInt = round(mean(interSample));
[pxx_meas, f_meas] = pwelch(irregSignal, hanning(avgInt), avgInt/2, [], targetSrate);
[pxx_up,   f_up]   = pwelch(newsignal,   hanning(winLen), winLen/2, [], targetSrate);

figure('Color','w','Position',[100 100 1000 450]); clf; hold on;
semilogy(f,      pxx,      'k-', 'LineWidth',2);
semilogy(f_up,   pxx_up,   'm-', 'LineWidth',1.5);
semilogy(f_meas, pxx_meas, 'r--','LineWidth',1.2);

% ---- MARKER AT THE *EXACT* BIN ----
plot(peakfreqBin, peakPower, 'kp','MarkerSize',14,'MarkerFaceColor','y',...
    'MarkerEdgeColor','k');

xlabel('Frequency (Hz)'); ylabel('Power (\muV^2/Hz)');
title(sprintf('PSD – Peak @ %.2f Hz (%s)', peakfreq, ...
    iif(useManualPeak,'Manual','Auto')));
legend({'Real','Reconstructed','Irregular'},'Location','best');
xlim([4 15]); grid on;

% Annotation
text(peakfreqBin, peakPower, sprintf(' %.2f Hz', peakfreqBin), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','bold','Color','k');

%% ----------------------- 8. Quantitative Metrics --------------------
R   = corr(signal, newsignal);
MAE = mean(abs(signal - newsignal));

fprintf('\n=== RESULTS ===\n');
fprintf('Peak Frequency (plotted): %.2f Hz (%s)\n', peakfreqBin, ...
    iif(useManualPeak,'Manual','Auto'));
fprintf('Correlation (real ↔ rec): %.4f\n', R);
fprintf('MAE: %.2f µV\n', MAE);

function out = iif(cond, t, f)
    if cond, out = t; else, out = f; end
end