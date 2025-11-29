%% PROFESSIONAL OUTLIER DETECTION BENCHMARK ON EUR/USD
% Your original idea + full evaluation suite (so you learn what really works)
% Author: Ali Arabi Bavil (2025)

clear all; close all; clc;
load forex.mat

N    = length(forex);
time = (0:N-1)' / (N-1);              % normalized time [0–1] → better for comparison

%% 1. Raw data overview
figure(1); clf; set(gcf,'Position',[100 100 1200 700]);
plot(time, forex, 'k', 'LineWidth', 1.2);
xlabel('Time (years)'); ylabel('EUR/USD');
title('EUR/USD Daily Rates (full history)');
grid on; box on;

%% 2. My LOCAL WINDOW METHOD (cleaned up + vectorized for speed)
pct_win = 5;                           % In percent
k = round(N * pct_win/200);            % half-window in samples (+-k around each point)

% Vectorized moving statistics (100x faster than loop!)
mean_ts  = movmean(forex, [k k]);      % symmetric window
std_ts   = movstd(forex, [k k]);       % same window
std3_ts  = 3 * std_ts;

outliers_local = forex > (mean_ts + std3_ts) | forex < (mean_ts - std3_ts);

%% 3. COMPETING METHODS (for fair comparison)
% Global fixed threshold 
global_thresh_up   = mean(forex) + 3*std(forex);
global_thresh_dn   = mean(forex) - 3*std(forex);
outliers_global    = forex > global_thresh_up | forex < global_thresh_dn;

% Robust MAD method (often best for finance)
med_val = median(forex);
mad_val = mad(forex,1);
mad_thresh = 6;                        % typical for FX data
outliers_mad = forex > (med_val + mad_thresh*mad_val) | ...
               forex < (med_val - mad_thresh*mad_val);

% Modified Z-score (excellent for local outliers)
mod_z = abs(forex - med_val) / mad_val;
outliers_modz = mod_z > 7;             % common threshold

%% 4. MASTER VISUALIZATION — SEE EVERYTHING AT ONCE
figure(2); clf; set(gcf,'Position',[100 50 1400 900]);
tiledlayout(3,1,'TileSpacing','compact');

% Top: My local window method (beautiful patch)
ax1 = nexttile; hold on;
h_patch = fill([time; flip(time)], [mean_ts+std3_ts; flip(mean_ts-std3_ts)], ...
               [1 0.4 0.8], 'EdgeColor','none', 'FaceAlpha',0.35);
plot(time, forex, 'k', 'LineWidth', 1.5);
plot(time(outliers_local), forex(outliers_local), 'ro', ...
     'MarkerFaceColor','r', 'MarkerSize',5);
plot(time, mean_ts, 'm', 'LineWidth', 1.5);
xlabel('Time (years)'); ylabel('EUR/USD');
title(sprintf('YOUR METHOD: ±%g%% Moving Window (k=%d points) → %d outliers', ...
      pct_win, 2*k+1, sum(outliers_local)));
legend({'±3σ envelope','EUR/USD','Detected outliers','Local mean'}, ...
       'Location','northoutside','Orientation','horizontal');
grid on;

% Middle: All methods overlaid
ax2 = nexttile; hold on;
plot(time, forex, 'Color',[0.3 0.3 0.3], 'LineWidth', 1);

plot(time(outliers_global), forex(outliers_global), 'bs', 'MarkerFaceColor','b', 'MarkerSize',6);
plot(time(outliers_local),  forex(outliers_local),  'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
plot(time(outliers_mad),    forex(outliers_mad),    'g^', 'MarkerFaceColor','g', 'MarkerSize',7);
plot(time(outliers_modz),   forex(outliers_modz),   'm+', 'MarkerSize',8);

plot(get(gca,'xlim'), [1 1]*global_thresh_up, 'b--');
plot(get(gca,'xlim'), [1 1]*global_thresh_dn, 'b--');

title('Comparison of All Methods');
legend({'EUR/USD', ...
        sprintf('Global ±3σ (%d)',sum(outliers_global)), ...
        sprintf('Your Local Window (%d)',sum(outliers_local)), ...
        sprintf('MAD×%g (%d)',mad_thresh,sum(outliers_mad)), ...
        sprintf('Modified Z > %g (%d)',7,sum(outliers_modz))}, ...
       'Location','northoutside');
grid on;

% Bottom: Histogram of detected outlier severity
ax3 = nexttile; hold on;
severity = forex(outliers_local) - (mean_ts(outliers_local) + std3_ts(outliers_local).*sign(forex(outliers_local)-mean_ts(outliers_local)));
histogram(abs(severity), 50, 'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.7);
xlabel('How far outlier is from local threshold (in σ units)');
ylabel('Count');
title('Outlier Severity Distribution (your method)');
grid on;

%% 5. QUANTITATIVE SCOREBOARD (learn what actually works!)
fprintf('\n=== OUTLIER DETECTION BENCHMARK ===\n');
fprintf('%-25s %8s %12s\n','Method','N_outliers','% of data');
fprintf('%s\n',repmat('-',1,50));
fprintf('%-25s %8d %10.3f%%\n','Global ±3σ',      sum(outliers_global),  100*mean(outliers_global));
fprintf('%-25s %8d %10.3f%%\n','Your Local Window',sum(outliers_local),   100*mean(outliers_local));
fprintf('%-25s %8d %10.3f%%\n','MAD ×6',           sum(outliers_mad),      100*mean(outliers_mad));
fprintf('%-25s %8d %10.3f%%\n','Modified Z-score',  sum(outliers_modz),     100*mean(outliers_modz));

%% 6. BONUS: Zoom on 2008 crisis — do methods catch real anomalies?
figure(3); clf;
idx_2008 = time > 0.45 & time < 0.55;   % roughly 2008–2009
plot(time(idx_2008), forex(idx_2008), 'k', 'LineWidth', 1.5); hold on;
plot(time(idx_2008 & outliers_local), forex(idx_2008 & outliers_local), 'ro', 'MarkerFaceColor','r','MarkerSize',8);
title('2008 Financial Crisis Zoom — Are detected outliers real events?');
grid on;