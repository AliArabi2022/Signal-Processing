%% EEG Upsampling: From 10 Hz Wearable to 100 Hz Clinical
% ======================================================
% - Real EEG from PhysioNet (downsampled to 10 Hz)
% - Upsample to 100 Hz using 4 methods
% - Compare vs. 1000 Hz ground truth
% - Metrics: MSE, PSNR, spectral correlation
clc; clear; close all;

%% 1. LOAD EDF
[hdr, events] = edfread('S001R03.edf');

%% 2. GET CHANNEL NAMES (EXCLUDE 'Record Time')
channel_names = string(hdr.Properties.VariableNames(2:end));  % ← 64 channels only
n_epochs = height(hdr);          % 125
samples_per_epoch = 160;
n_channels = numel(channel_names);  % ← 64, NOT width(hdr)!

%% 3. PREALLOCATE
total_samples = n_epochs * samples_per_epoch;
eeg_data = zeros(total_samples, n_channels);

%% 4. EXTRACT DATA
for ch = 1:n_channels
    col_name = channel_names(ch);           % ← string, safe
    epoch_data = hdr.(col_name);            % ← 125×1 cell
    for ep = 1:n_epochs
        start_idx = (ep-1)*samples_per_epoch + 1;
        end_idx = ep * samples_per_epoch;   % ← defined!
        eeg_data(start_idx:end_idx, ch) = epoch_data{ep};
    end
end

%% 5. TIME VECTOR
fs = 160;
t_full = (0:total_samples-1) / fs;

%% 6. EVENTS
event_times = seconds(events.Onset);
event_labels = events.Annotations;

%% 7. PICK Cz
cz_idx = find(channel_names == "Cz__");
if isempty(cz_idx)
    error('Cz__ not found! Channels: %s', join(channel_names, ', '));
end
eeg_cz = eeg_data(:, cz_idx);

%% 8. PLOT Cz + EVENTS
figure('Position', [100, 100, 1200, 500]);
plot(t_full, eeg_cz, 'Color', [0.1 0.4 0.7], 'LineWidth', 1.1); hold on;

colors = containers.Map({'T0','T1','T2'}, {'g','r','b'});
for i = 1:length(event_labels)
    t_evt = event_times(i);
    label = event_labels{i};
    if isKey(colors, label)
        idx = round(t_evt * fs) + 1;
        if idx >= 1 && idx <= length(eeg_cz)
            plot(t_evt, eeg_cz(idx), 'o', ...
                'MarkerFaceColor', colors(label), 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
            text(t_evt, max(eeg_cz)*1.15, label, ...
                'Color', colors(label), 'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'FontSize', 10);
        end
    end
end

xlabel('Time (s)'); ylabel('Amplitude (µV)');
title('Cz Channel – Real 64-Ch EEG (160 Hz) with Motor Events');
grid on; axis tight;
legend('EEG (Cz)', 'Events');


%% 9. SEGMENT AROUND FIRST T1
t1_idx = find(strcmp(event_labels, 'T1'), 1);
seg_start = event_times(t1_idx) - 1;
seg_end = seg_start + 5;

idx_seg = t_full >= seg_start & t_full <= seg_end;
eeg_seg = eeg_cz(idx_seg);
t_seg = t_full(idx_seg);

% Downsample to 10 Hz
fs_low = 10;
dec_factor = fs / fs_low;
idx_low = 1:dec_factor:length(eeg_seg);
eeg_low = eeg_seg(idx_low);
t_low = t_seg(idx_low);

% Target 100 Hz
fs_target = 100;
t_target = linspace(t_low(1), t_low(end), round(5 * fs_target) + 1);

gt_100 = interp1(t_seg, eeg_seg, t_target, 'spline');
eeg_up = interp1(t_low, eeg_low, t_target, 'spline');

%% 10. PLOT
figure;
plot(t_low, eeg_low, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8); hold on;
plot(t_target, gt_100, 'k--', 'LineWidth', 1.8);
plot(t_target, eeg_up, 'r-', 'LineWidth', 1.4);
legend('10 Hz', '100 Hz GT', 'Upsampled');
title('Real EEG: 10 Hz → 100 Hz (T1 Motor Task)');
grid on;