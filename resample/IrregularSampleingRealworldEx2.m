
clear; close all; clc;

%% ----------------------- 1. Simulation parameters --------------------
srate    = 1324;               % Hz
duration = 2;                  % seconds
npnts    = srate * duration;   % total samples
timevec  = (0:npnts-1)/srate;  % regular time vector (s)
peakfreq = 12;                  % Hz – centre of alpha
fwhm     = 5;                  % Hz – full-width at half-maximum
hz       = linspace(0,srate,npnts);   % frequency vector for FFT

%% ----------------------- 2. Build the clean signal -------------------
s = fwhm*(2*pi-1)/(4*pi);               % normalized width
x = hz - peakfreq;                      % shifted frequencies
fg = exp(-.5*(x/s).^2);                 % Gaussian taper

fc = rand(1,npnts) .* exp(1i*2*pi*rand(1,npnts));
fc = fc .* fg;                          % apply alpha peak

signal = 2 * real(ifft(fc)) * npnts;    % inverse FFT → time domain
signal = signal(:);                     % column vector

%% ----------------------- 3. Irregular sampling ----------------------
rng(42);                                % reproducibility
interSample = ceil(exp(4*rand(npnts,1)));
sampIdx = cumsum([1; interSample]);
sampIdx(sampIdx > npnts) = [];          % clip overflow

irregSignal = signal(sampIdx);
irregTime   = timevec(sampIdx);

% ---- SORT BY TIME (required for griddedInterpolant) ----
[irregTime, sortIdx] = sort(irregTime);
irregSignal = irregSignal(sortIdx);

%% ----------------------- 4. Upsample to original grid ---------------
F = griddedInterpolant(irregTime, irregSignal, 'spline');
newsignal = F(timevec);                 % evaluate on regular grid
newsignal = newsignal(:);

%% ----------------------- 5. Visualisation (time domain) ----------
figure('Color','w','Position',[100 100 900 600]); clf; hold on;
plot(timevec, signal,      'k', 'LineWidth',2);               % Analog
plot(irregTime, irregSignal, 'ro', 'MarkerFaceColor','r','MarkerSize',8); % Measured
plot(timevec, newsignal,   'm.', 'MarkerSize',6);            % Upsampled

xlabel('Time (s)'); ylabel('Amplitude (a.u.)');
title('Irregularly sampled EEG → reconstruction');
legend({'"Analog" (1324 Hz)','Measured (irregular)','Upsampled (spline)'},...
       'Location','southoutside');
grid on; axis tight;

%% ----------------------- 6. Spectral comparison --------------------
winLen = srate;                         % 1-second window
noverlap = winLen / 2;

% Analog
[pxx_analog, f] = pwelch(signal, hanning(winLen), noverlap, [], srate);

% Irregular: use average inter-sample distance to define window
avgInterval = round(mean(interSample));
winIrreg = hanning(avgInterval);
noverlapIrreg = round(avgInterval / 2);

[pxx_meas, f_meas] = pwelch(irregSignal, winIrreg, noverlapIrreg, [], srate);

% Upsampled
[pxx_up, f_up] = pwelch(newsignal, hanning(winLen), noverlap, [], srate);

% Plot
figure('Color','w','Position',[100 100 900 400]); clf; hold on;
semilogy(f,      pxx_analog, 'k',  'LineWidth',2);
semilogy(f_up,   pxx_up,     'm',  'LineWidth',1.5);
semilogy(f_meas, pxx_meas,   'r--','LineWidth',1.2);
xlabel('Frequency (Hz)'); ylabel('Power (a.u.^2/Hz)');
title('Power Spectral Density');
legend({'Analog','Upsampled','Measured (irregular)'}, 'Location','best');
xlim([0 30]); grid on;
%% ----------------------- 7. Quantitative check ----------------------
R   = corr(signal, newsignal);
MAE = mean(abs(signal - newsignal));

fprintf('\nPearson correlation (analog ↔ upsampled): %.4f\n', R);
fprintf('Mean absolute error: %.6f a.u.\n', MAE);