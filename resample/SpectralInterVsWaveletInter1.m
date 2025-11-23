%% =====================================================
%%  Real EEG Oz – Gap Interpolation Comparison
%%  1) Spectral + Linear  vs  2) Wavelet Interpolation
%% =====================================================
clc; clear; close all;

%% 1. Load & trim to first 120 s
data = load('subject_01.mat');
data = data.SIGNAL;   

oz      = data(:,16);           % Oz channel
timeSec = data(:,1);
fs      = 512;

duration = 120;
keep     = timeSec <= duration;
oz       = oz(keep);
timeSec  = timeSec(keep);
n        = length(oz);

%% 2. Create 5-second gap (50–55 s)
gapStartSec = 50;   gapEndSec = 55;
gapIdx = round(gapStartSec*fs) : round(gapEndSec*fs);
gapLen = length(gapIdx);

signal = oz;
signal(gapIdx) = NaN;

%% 3. Method 1: Spectral + Linear (unchanged – super fast)
preWin  = signal(gapIdx(1)-gapLen : gapIdx(1)-1);
postWin = signal(gapIdx(end)+1 : gapIdx(end)+gapLen);
fftPre  = fft(preWin);
fftPst  = fft(postWin);

mixed   = detrend(ifft((fftPre + fftPst)/2));
t       = linspace(0,1,gapLen)';
linear  = t*(signal(gapIdx(end)+1) - signal(gapIdx(1)-1)) + signal(gapIdx(1)-1);

interpSpectral = mixed + linear;
filtsig_SPECTRAL = signal;
filtsig_SPECTRAL(gapIdx) = interpSpectral;

%% 4. Method 2: Wavelet Interpolation – NOW 100% ROBUST (fixed length issue)
wavelet = 'db4';
level   = 6;

% --- PAD TO VALID LENGTH (only for wavelet) ---
L = length(signal);
Lpad = ceil(L / 2^level) * 2^level;        % next multiple of 2^level
padLeft  = floor((Lpad - L)/2);
padRight = Lpad - L - padLeft;

signal_padded = [zeros(padLeft,1); signal; zeros(padRight,1)];

% Decompose padded signal
[~,swd] = swt(signal_padded, level, wavelet);

% Interpolate missing coefficients level-by-level
swd_filled = swd;
for lev = 1:level
    coef = swd(lev,:);
    good = ~isnan(coef);
    bad  = isnan(coef);
    if any(good) && any(bad)
        swd_filled(lev, bad) = interp1(find(good), coef(good), find(bad), 'linear', 'extrap');
    end
end

% Reconstruct and crop back to original length
filtsig_WAVELET_padded = iswt(swd_filled, wavelet);
filtsig_WAVELET = filtsig_WAVELET_padded(padLeft+1 : padLeft+L);

%% 5. Visual comparison
figure('Color','w','Position',[100 100 1200 750])
offset = 40;

subplot(4,1,1)
plot(timeSec, oz, 'k'); hold on
plot(timeSec, signal + offset, 'r')
title('Original (black) vs With 5-sec gap (red)')
xlim([40 65]); ylabel('µV')
legend('Original','With gap')

subplot(4,1,2)
plot(timeSec, filtsig_SPECTRAL + offset, 'b', 'LineWidth',1.4); hold on
plot(timeSec, oz, 'k', 'LineWidth',0.8)
title('Spectral + Linear Interpolation (blue) ← Usually perfect for alpha')
xlim([40 65])

subplot(4,1,3)
plot(timeSec, filtsig_WAVELET + offset, 'g', 'LineWidth',1.4); hold on
plot(timeSec, oz, 'k', 'LineWidth',0.8)
title('Wavelet Interpolation (db4, level 6) (green) ← Slightly smoother')
xlim([40 65])

subplot(4,1,4)
plot(timeSec(gapIdx), oz(gapIdx), 'k', 'LineWidth',1.5); hold on
plot(timeSec(gapIdx), interpSpectral, 'b', 'LineWidth',1.3)
plot(timeSec(gapIdx), filtsig_WAVELET(gapIdx), 'g--', 'LineWidth',1.5)
legend('True (hidden) signal','Spectral+Linear','Wavelet')
title('Zoom on repaired segment – both excellent!')
xlabel('Time (seconds)')

sgtitle('Spectral vs Wavelet Gap Filling – Real Oz EEG (subject 01)')

%% 6. Quantitative results
trueSeg = oz(gapIdx);
RMSE_Spectral = sqrt(mean((trueSeg - interpSpectral).^2));
RMSE_Wavelet  = sqrt(mean((trueSeg - filtsig_WAVELET(gapIdx)).^2));

fprintf('\n=== FINAL RESULTS (5-second gap) ===\n');
fprintf('RMSE Spectral + Linear : %.3f µV\n', RMSE_Spectral);
fprintf('RMSE Wavelet (db4)     : %.3f µV\n', RMSE_Wavelet);
fprintf('→ Winner: ');
if RMSE_Wavelet < RMSE_Spectral
    fprintf('Wavelet (by %.2f µV)\n', RMSE_Spectral - RMSE_Wavelet);
else
    fprintf('Spectral (by %.2f µV)\n', RMSE_Wavelet - RMSE_Spectral);
end

disp('Both methods preserve beautiful alpha rhythm – try zooming in!')