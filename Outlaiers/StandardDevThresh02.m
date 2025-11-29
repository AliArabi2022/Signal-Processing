%% Real-World Outlier Detection & Cleaning on Red Wine Quality Dataset
% Dataset: https://archive.ics.uci.edu/ml/datasets/wine+quality
% File:    winequality-red.csv
%
% We will treat one of the most outlier-prone columns ("residual sugar" or
% "chlorides" or "sulphates") as our time-series-like signal and clean it
% using the exact same robust MAD + bidirectional interpolation method.
%
% Author: Ali Arabi Bavil
% Date:   2025

clear; close all; clc;

%% 1. Load the real data
data = readtable('winequality-red.csv');   % Make sure file is in current folder
disp(data(1:6,:));                         % Peek at first few rows

%% 2. Choose a column that typically has outliers
% These three are the most common ones with heavy tails / spikes:
%   - residual sugar   → some wines are extremely sweet
%   - chlorides        → salt content, occasional very high values
%   - sulphates        → preservative, some extreme values in cheap wines

% signal_raw = data.sulphates;           % ← TRY ALSO: data.residualSugar or data.chlorides
% variableName = 'sulphates';            % Change accordingly if you switch

% signal_raw = data.residualSugar;
% variableName = 'residual sugar';  % very heavy tail!

signal_raw = data.chlorides;
variableName = 'chlorides';        % many spikes

% signal_raw = data.totalSulfurDioxide;
% variableName = 'total sulfur dioxide';



N = length(signal_raw);
time = (1:N)';                         % Treat samples as sequential (like a time series)

fprintf('\nAnalyzing column: %s (%d samples)\n', variableName, N);

%% 3. Plot raw data with obvious outliers
figure(1); clf; hold on;
plot(time, signal_raw, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.6 0.6 0.6], ...
     'DisplayName', 'Raw data');

xlabel('Sample Index (Wine #)');
ylabel(sprintf('%s (g/dm³)', variableName));
title(sprintf('Red Wine Dataset – Raw %s (with outliers)', upper(variableName)));
grid on;

%% 4. Robust bidirectional outlier detection using MAD (perfect for skewed data)

medVal = median(signal_raw);
madVal = mad(signal_raw, 1);           % 1 → already scaled for normal distribution

k = 6;                                 % 5–7 works great on this dataset (tune if needed)

lowerThresh = medVal - k * madVal;
upperThresh = medVal + k * madVal;

% Plot thresholds
plot([1 N], [upperThresh upperThresh], 'r--', 'LineWidth', 2, 'DisplayName', 'Upper threshold');
plot([1 N], [lowerThresh lowerThresh], 'r--', 'LineWidth', 2, 'DisplayName', 'Lower threshold');

fprintf('Median = %.4f   |   MAD = %.4f\n', medVal, madVal);
fprintf('Thresholds:  lower = %.4f    upper = %.4f\n', lowerThresh, upperThresh);

%% 5. Find outliers (both sides)
outliers = signal_raw < lowerThresh | signal_raw > upperThresh;
nOutliers = sum(outliers);

fprintf('%d outliers detected (%.2f%% of total samples)\n', nOutliers, 100*nOutliers/N);

% Highlight them
plot(time(outliers), signal_raw(outliers), 'rp', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Detected outliers');

%% 6. Clean by linear interpolation (exactly like your original code)
goodPoints = ~outliers;

F = griddedInterpolant(time(goodPoints), signal_raw(goodPoints), ...
                       'linear', 'nearest');   % 'nearest' prevents NaN at edges

signal_clean = signal_raw;
signal_clean(outliers) = F(time(outliers));

% Plot cleaned version
plot(time, signal_clean, 'b-', 'LineWidth', 2, 'DisplayName', 'Cleaned signal');

legend('Location', 'best');
title(sprintf('%s – Outliers Removed via Interpolation (%d replaced)', ...
      upper(variableName), nOutliers));

%% 7. Before / After comparison (zoomed view of worst region
figure(2); clf;

subplot(2,1,1);
plot(time, signal_raw, 'k.-', 'MarkerSize', 4);
hold on;
plot(time(outliers), signal_raw(outliers), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot([1 N], [upperThresh upperThresh], 'r--');
plot([1 N], [lowerThresh lowerThresh], 'r--');
ylabel(sprintf('%s (raw)', variableName));
title('Before Cleaning');
grid on;

subplot(2,1,2);
plot(time, signal_clean, 'b.-', 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(sprintf('%s (cleaned)', variableName));
xlabel('Sample Index (Wine #)');
title(sprintf('After Cleaning – %d outliers interpolated', nOutliers));
grid on;

sgtitle(sprintf('Wine Quality Dataset – %s Outlier Removal', upper(variableName)));

%% 8. Statistics summary
disp(' ');
disp('=== Cleaning Summary ===');
disp(table([medVal; median(signal_clean)], ...
           [madVal; mad(signal_clean,1)], ...
           [max(signal_raw); max(signal_clean)], ...
           [min(signal_raw); min(signal_clean)], ...
    'VariableNames', {'Median', 'MAD', 'Max', 'Min'}, ...
    'RowNames', {'Raw', 'Cleaned'}));

% Optional: save cleaned dataset
data.([variableName '_clean']) = signal_clean;
% writetable(data, 'winequality-red_cleaned.csv');   % Uncomment to save