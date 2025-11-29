%% Real-World Outlier Detection & Cleaning – Red Wine Quality Dataset
% Now with side-by-side histograms (before vs after cleaning)

clear; close all; clc;

%% 1. Load data
data = readtable('winequality-red.csv');

%% 2. Choose column (most interesting ones for outliers)
signal_raw   = data.sulphates;           % ← try also: residualSugar, chlorides, totalSulfurDioxide
variableName = 'sulphates';
unit         = 'g/dm³';

N = length(signal_raw);
time = (1:N)';

fprintf('Analyzing: %s (%d samples)\n\n', variableName, N);

%% 3. Robust MAD-based bidirectional outlier detection
medVal = median(signal_raw);
madVal = mad(signal_raw, 1);

k = 6;                                        % 5–7 works perfectly on this dataset
lowerThresh = medVal - k * madVal;
upperThresh = medVal + k * madVal;

outliers = signal_raw < lowerThresh | signal_raw > upperThresh;
nOutliers = sum(outliers);

fprintf('Thresholds:  %.4f  to  %.4f\n', lowerThresh, upperThresh);
fprintf('%d outliers detected (%.2f%%)\n\n', nOutliers, 100*nOutliers/N);

%% 4. Clean via linear interpolation
F = griddedInterpolant(time(~outliers), signal_raw(~outliers), 'linear', 'nearest');

signal_clean = signal_raw;
signal_clean(outliers) = F(time(outliers));

%% 5. Plot time-series (before & after)
figure(1); clf;
subplot(4,1,[1 2]); hold on; box on;

plot(time, signal_raw, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', [0.7 0.7 0.7]);
plot(time(outliers), signal_raw(outliers), 'rp', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(time, signal_clean, 'b-', 'LineWidth', 1.8);

plot([1 N], [upperThresh upperThresh], 'r--', 'LineWidth', 1.2);
plot([1 N], [lowerThresh lowerThresh], 'r--', 'LineWidth', 1.2);

xlabel('Sample Index (Wine #)');
ylabel({upper(variableName), ['(' unit ')']});
title(sprintf('%s – Outlier Removal via Interpolation (%d outliers replaced)', ...
    upper(variableName), nOutliers));
legend({'Raw data', 'Detected outliers', 'Cleaned signal', 'Thresholds'}, ...
       'Location', 'best');
grid on;

%% 6. HISTOGRAM COMPARISON – Raw vs Cleaned
subplot(4,1,[3,4]); hold on; box on;

% Use same bins for fair comparison
edges = linspace(min(signal_raw), max(signal_raw), 50);

histogram(signal_raw, edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', ...
          'FaceAlpha', 0.6, 'DisplayName', 'Raw data');
histogram(signal_clean, edges, 'FaceColor', [0 0.45 0.74], 'EdgeColor', 'k', ...
          'FaceAlpha', 0.7, 'DisplayName', 'Cleaned data');

% Overlay detected outliers as red bars on top
if nOutliers > 0
    histogram(signal_raw(outliers), edges, 'FaceColor', 'r', 'EdgeColor', 'k', ...
              'FaceAlpha', 0.8, 'DisplayName', 'Outliers');
end

xlabel({upper(variableName), ['(' unit ')']});
ylabel('Count');
title('Histogram: Before vs After Outlier Removal');
legend('Location', 'northoutside', 'Orientation', 'horizontal');
grid on;

% Make the whole figure nice and tall
set(gcf, 'Position', [100 100 900 750]);

%% 7. Bonus: Inset zoomed histogram of the tail (optional – uncomment if you like)
axes('Position',[0.55 0.22 0.35 0.25]);
histogram(signal_raw, 100, 'FaceColor', [0.8 0.8 0.8], 'DisplayName', 'Raw');
histogram(signal_clean, 100, 'FaceColor', [0 0.4 0.8], 'DisplayName', 'Cleaned');
xlim([upperThresh-0.2 max(signal_raw)]);
title('Zoom on upper tail');
grid on;

%% 8. Summary table
disp('=== Cleaning Summary ===');
disp(table([medVal; median(signal_clean)], ...
           [madVal; mad(signal_clean,1)], ...
           [max(signal_raw); max(signal_clean)], ...
           [min(signal_raw); min(signal_clean)], ...
           [nOutliers; 0], ...
    'VariableNames', {'Median', 'MAD', 'Max', 'Min', 'N_Outliers'}, ...
    'RowNames', {'Raw', 'Cleaned'}));

% Optional: save cleaned version
data.([variableName '_clean']) = signal_clean;
% writetable(data, 'winequality-red_with_clean_sulphates.csv');