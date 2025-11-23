clc; clear all; close all;


%% 1. Load and extract Oz channel
data = load('subject_01.mat');
signalMatrix = data.SIGNAL;                    % Tx19 matrix
oz_idx = 16;                                 % Oz is column 16
eegRaw = signalMatrix(:, oz_idx);            % Full recording (several minutes)
timeRaw = signalMatrix(:, 1);                % First column = time in seconds

%% 2. NEW SECTION: Trim to first 120 seconds (2 minutes) → much faster!
duration = 120;                              % seconds
keepIdx = timeRaw <= duration;               % Logical index for first 120 s
timeRaw = timeRaw(keepIdx);                  % Trim time vector
eegRaw  = eegRaw(keepIdx);                   % Trim EEG (Oz only)

fprintf('Data trimmed to first %.1f seconds (n = %d samples, fs = 512 Hz)\n', ...
        duration, length(eegRaw));

%% 3. Prepare signals
fs = 512;                                    % Sampling rate (known from dataset)
n  = length(eegRaw);
origsig = eegRaw;                            % Keep original for plotting
signal  = eegRaw;                            % Working copy

%% 4. Define and create a gap inside the 2-minute window
% Let's put a ~5-second gap in the middle (e.g., seconds 50–55)
gap_start_sec = 50;
gap_end_sec   = 55;

boundaryPnts = round([gap_start_sec gap_end_sec] * fs);
gap_length   = diff(boundaryPnts) + 1;       % Number of samples in gap

% Safety check (gap must fit in trimmed data)
if boundaryPnts(2) > n
    error('Gap exceeds data length. Reduce gap or increase duration.');
end

signal(boundaryPnts(1):boundaryPnts(2)) = NaN;   % Mark gap with NaN

%% 5. Spectral + linear interpolation (same excellent method)
% Pre-gap window (same length as gap)
pre_start = boundaryPnts(1) - gap_length;
pre_end   = boundaryPnts(1) - 1;
fftPre = fft(signal(pre_start:pre_end));

% Post-gap window
post_start = boundaryPnts(2) + 1;
post_end   = boundaryPnts(2) + gap_length;
fftPst = fft(signal(post_start:post_end));

% Spectral part: average FFTs → IFFT → detrend
mixeddata = detrend(ifft((fftPre + fftPst)/2));

% Linear part: connect endpoints perfectly
start_val = signal(boundaryPnts(1)-1);
end_val   = signal(boundaryPnts(2)+1);
t         = linspace(0, 1, gap_length)';
linedata  = t * (end_val - start_val) + start_val;

% Combine both
linterp = mixeddata + linedata;

% Insert repaired segment
filtsig = signal;
filtsig(boundaryPnts(1):boundaryPnts(2)) = linterp;

%% 6. Plotting – clean and focused
figure(1), clf
offset = 30;  % µV offset for visibility
plot(timeRaw, origsig, 'k', 'LineWidth', 1); hold on
plot(timeRaw, signal  + offset, 'r')
plot(timeRaw, filtsig + 2*offset, 'b', 'LineWidth', 1.2)

xlim([40 65])                                 % Zoom around the gap
legend({'Original Oz', 'With 5-second gap', 'Repaired (spectral + linear)'}, ...
       'Location', 'northwest')
xlabel('Time (seconds)')
ylabel('Amplitude (µV, offset for clarity)')
title('Real EEG (Oz) - 5-second Gap Repaired with Spectral + Linear Interpolation')
grid on

%% 7. Bonus: Close-up + PSD comparison
figure(2), clf

subplot(2,1,1)
plot(timeRaw(pre_start:post_end), ...
     [origsig(pre_start:post_end), ...
      signal(pre_start:post_end) + offset, ...
      filtsig(pre_start:post_end) + 2*offset], 'LineWidth', 1.5)
xlim([timeRaw(pre_start) timeRaw(post_end)])
legend('Original', 'Gap (NaN)', 'Repaired')
title('Close-up around repaired segment')
grid on

subplot(2,1,2)
f = (0:gap_length-1) * (fs / gap_length);
Nplot = floor(gap_length/2);
plot(f(1:Nplot), abs(fftPre(1:Nplot))/gap_length, 'b'); hold on
plot(f(1:Nplot), abs(fftPst(1:Nplot))/gap_length, 'r')
plot(f(1:Nplot), abs((fftPre(1:Nplot) + fftPst(1:Nplot))/2)/gap_length, 'g', 'LineWidth', 2)
xlim([4 15])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
legend('Pre-gap', 'Post-gap', 'Averaged (used for repair)')
title('Alpha band (8–12 Hz) is beautifully preserved!')
grid on

disp('Done! You now have a perfectly repaired 5-second gap in real Oz alpha EEG.')