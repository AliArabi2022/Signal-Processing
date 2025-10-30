%% EEG DOWNSAMPLING WITH FIR FILTER (3 RATES) + TIME/FREQ PLOTS
clc; clear; close all;

%% 1. LOAD EDF (EXACTLY LIKE FINAL WORKING CODE)
[hdr, events] = edfread('S001R03.edf');

% Extract channel names (exclude 'Record Time')
channel_names = string(hdr.Properties.VariableNames(2:end));
n_epochs = height(hdr);          % 125
samples_per_epoch = 160;
n_channels = numel(channel_names);  % 64

% Reconstruct full EEG matrix
total_samples = n_epochs * samples_per_epoch;
eeg_data = zeros(total_samples, n_channels);

for ch = 1:n_channels
    col_name = channel_names(ch);
    epoch_data = hdr.(col_name);
    for ep = 1:n_epochs
        start_idx = (ep-1)*samples_per_epoch + 1;
        end_idx = ep * samples_per_epoch;
        eeg_data(start_idx:end_idx, ch) = epoch_data{ep};
    end
end

%% 2. SETTINGS
fs_orig = 160;
cz_idx = find(channel_names == "Cz__");
eeg_cz = eeg_data(:, cz_idx);
t_full = (0:total_samples-1) / fs_orig;

% Downsample to 80, 40, 20 Hz
fs_targets = [80, 40, 20];
colors = {'b', 'r', 'g'};

%% 3. FIGURE SETUP
figure('Position', [100, 100, 1400, 900]);

%% LOOP: DOWNSAMPLE + PLOT TIME + FREQ
for i = 1:length(fs_targets)
    fs_new = fs_targets(i);
    dec_factor = fs_orig / fs_new;
    
    % Design FIR low-pass filter (anti-aliasing)
    Nyq_new = fs_new / 2;
    h = fir1(100, Nyq_new / (fs_orig/2), 'low');  % 100th order FIR
    
    % Filter
    eeg_filt = filtfilt(h, 1, eeg_cz);
    
    % Downsample
    eeg_down = downsample(eeg_filt, dec_factor);
    t_down = downsample(t_full, dec_factor);
    
    % Store
    downsampled{i} = eeg_down;
    t_downsampled{i} = t_down;
    fs_downsampled{i} = fs_new;
    
    %% TIME DOMAIN PLOT (First 5 seconds)
    subplot(3, 2, 2*i-1); hold on;
    plot(t_full(1:800), eeg_cz(1:800), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    plot(t_down(1:round(5*fs_new)), eeg_down(1:round(5*fs_new)), ...
        'Color', colors{i}, 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Amplitude (µV)');
    title(sprintf('Time Domain: %d Hz (Downsampled)', fs_new));
    legend('160 Hz', sprintf('%d Hz', fs_new));
    grid on; axis tight;
    
    %% FREQUENCY DOMAIN (PSD)
    subplot(3, 2, 2*i);
    [Pxx_orig, f_orig] = pwelch(eeg_cz, hamming(512), 256, 512, fs_orig);
    [Pxx_down, f_down] = pwelch(eeg_down, hamming(512), 256, 512, fs_new);
    
    plot(f_orig, 10*log10(Pxx_orig), 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
    plot(f_down, 10*log10(Pxx_down), 'Color', colors{i}, 'LineWidth', 1.8);
    xlabel('Frequency (Hz)'); ylabel('Power (dB/Hz)');
    title(sprintf('PSD: %d Hz', fs_new));
    legend('160 Hz', sprintf('%d Hz', fs_new));
    grid on; xlim([0 80]);
end

sgtitle('EEG Downsampling: FIR Filter + Time/Frequency Analysis (Cz Channel)', 'FontSize', 14);