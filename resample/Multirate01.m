%% MULTIRATE EEG: 10, 40, 83 Hz → 160 Hz (CLEAR LEGENDS)
clc; clear; close all;

%% 1. LOAD S001R03.edf (YOUR WORKING METHOD)
[hdr, events] = edfread('S001R03.edf');
channel_names = string(hdr.Properties.VariableNames(2:end));
n_epochs = height(hdr); samples_per_epoch = 160;
n_channels = numel(channel_names);
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
t_full = (0:total_samples-1) / fs_orig;

% Pick 3 channels
chan_labels = ["Cz__", "C3__", "Fz__"];
chan_idx = arrayfun(@(x) find(channel_names == x), chan_labels);
eeg_selected = eeg_data(:, chan_idx);

% Target rates
fs_target = [10, 40, 83];
colors = 'kbr'; shapes = 'os^';

%% 3. DOWNSAMPLE USING resample()
[fs_down, time_down, signals_down] = deal(cell(3,1));
for si = 1:3
    [P, Q] = rat(fs_target(si) / fs_orig, 1e-6);
    sig_down = resample(eeg_selected(:, si), P, Q);
    t_down = (0:length(sig_down)-1) / fs_target(si);
    
    fs_down{si} = fs_target(si);
    time_down{si} = t_down;
    signals_down{si} = sig_down;
end

%% 4. FIGURE 1: DOWNSAMPLED (CLEAR LEGEND)
figure(1), clf, hold on;
legend_entries = {};
for si = 1:3
    h = plot(time_down{si}, signals_down{si}, ...
        [colors(si) shapes(si) '-'], 'LineWidth', 1.2, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6);
    legend_entries{end+1} = sprintf('%s (%.0f Hz)', chan_labels(si), fs_target(si));
end
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title('Downsampled EEG Signals');
grid on; axis tight;
legend(legend_entries, 'Location', 'northeastoutside', 'FontSize', 10);

%% 5. UPSAMPLE TO 160 Hz
newTime = t_full;
newsignals = zeros(3, length(newTime));
for si = 1:3
    F = griddedInterpolant(time_down{si}, signals_down{si}, 'spline');
    newsignals(si, :) = F(newTime);
end

%% 6. FIGURE 2: UPSAMPLED vs ORIGINAL (CLEAR LEGEND)
figure(2), clf, hold on;
legend_entries = {};
h_orig = [];
h_up = [];

for si = 1:3
    % Upsampled
    h_up(si) = plot(newTime, newsignals(si,:), [colors(si) '-'], ...
        'LineWidth', 2.2);
    legend_entries{end+1} = sprintf('%s: %.0f→160 Hz (Upsampled)', chan_labels(si), fs_target(si));
    
    % Original (gray background)
    h_orig(si) = plot(t_full, eeg_selected(:, si), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
end

xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title('Multirate Upsampling: Reconstruction to 160 Hz');
grid on; axis tight;

% Legend: Group upsampled + original
legend([h_up, h_orig(1)], ...
    [legend_entries, 'Original 160 Hz (Reference)'], ...
    'Location', 'northeastoutside', 'FontSize', 10);

%% 7. FIGURE 3: PSD (CLEAR LEGEND)
figure(3), clf;
for si = 1:3
    subplot(3,1,si); hold on;
    
    [Pxx_orig, f_orig] = pwelch(eeg_selected(:, si), [], [], 1024, fs_orig);
    [Pxx_up, f_up] = pwelch(newsignals(si,:), [], [], 1024, fs_orig);
    
    plot(f_orig, 10*log10(Pxx_orig), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
    plot(f_up, 10*log10(Pxx_up), 'Color', colors(si), 'LineWidth', 2.2);
    
    xlabel('Frequency (Hz)'); ylabel('Power (dB/Hz)');
    title(sprintf('%s Channel', chan_labels(si)));
    legend('Original 160 Hz', 'Upsampled', 'Location', 'northeast');
    grid on; xlim([0 80]);
end
sgtitle('Power Spectral Density: Original vs Upsampled', 'FontSize', 14, 'FontWeight', 'bold');