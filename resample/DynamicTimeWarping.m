%% Dynamic Time Warping (DTW) Demonstration
% Professional implementation comparing manual DTW vs MATLAB's built-in dtw()
% Author: Ali Arabi Bavil | Date: 2025
% Purpose: Educational + Real-world prototype for time series alignment
%
% Real-world applications:
%  • Speech recognition (word alignment under variable speed)
%  • Gesture recognition (smartphones, wearables, VR)
%  • Biomedical signal analysis (ECG, EEG, gait)
%  • Audio fingerprinting and music retrieval
%  • Anomaly detection in sensor data with timing variations

clear; clc; close all;

%% 1. Test Signals
tx = linspace(0, 1.5*pi, 400);   % Reference time vector
ty = linspace(0, 8*pi, 100);     % Query time vector

x = sin(tx.^2);                  % Chirp (accelerating frequency)
y = sin(ty);                     % Regular sine

%% 2. Original Signals
figure(1); clf;
subplot(3,1,1); hold on; grid on;
plot(tx, x, 'b.-', 'LineWidth',1.2);
plot(ty, y, 'rs-', 'LineWidth',1.4, 'MarkerFaceColor','w');
xlabel('Time (rad)'); ylabel('Amplitude');
title('Original Signals');
legend('Chirp x(t)', 'Sine y(t)');

%% 3. Manual DTW - Accumulated Cost Matrix
N = length(x);
M = length(y);

D = zeros(N, M);                                    % Accumulated cost
D(1,1) = abs(x(1) - y(1));

% First row
for i = 2:N
    D(i,1) = D(i-1,1) + abs(x(i) - y(1));
end

% First column
for j = 2:M
    D(1,j) = D(1,j-1) + abs(x(1) - y(j));
end

% Main DP fill
for i = 2:N
    for j = 2:M
        cost = abs(x(i) - y(j));                    % Try (x(i)-y(j))^2 for squared error
        D(i,j) = cost + min([ D(i-1,j), ...         % from above
                              D(i,j-1), ...         % from left
                              D(i-1,j-1) ]);        % from diagonal
    end
end

%% 4. SAFE BACKTRACKING - This version NEVER accesses invalid indices
i = N;
j = M;

% Pre-allocate with plenty of room
p = zeros(N+M, 1);   % indices in x
q = zeros(N+M, 1);   % indices in y
k = 1;

p(k) = i;
q(k) = j;
k = k + 1;

while i > 1 || j > 1                               % Stop only when we reach (1,1)
    if i == 1                                      % We are on the first row - can only go left
        j = j - 1;
    elseif j == 1                                  % We are on the first column - can only go up
        i = i - 1;
    else
        % Normal case: choose the best predecessor
        [~, idx] = min([D(i-1,j), D(i,j-1), D(i-1,j-1)]);
        switch idx
            case 1
                i = i - 1;
            case 2
                j = j - 1;
            case 3
                i = i - 1;
                j = j - 1;
        end
    end
    
    p(k) = i;
    q(k) = j;
    k = k + 1;
end

% Trim vectors to actual size and reverse (so we go from 1->N)
p = p(1:k-1);
q = q(1:k-1);
p = flipud(p);
q = flipud(q);

%% 5. Plot Manual DTW Result
subplot(3,1,2); hold on; grid on;
plot(tx, x, 'b-', 'LineWidth', 1.5);
plot(tx(p), y(q), 'r-', 'LineWidth', 2.2);
scatter(tx(p), y(q), 50, 'ro', 'filled');
xlabel('Time (aligned to x)');
title('Manual DTW - Perfect & Safe Alignment');
legend('Reference x(t)', 'Warped y(t) -> x time base');

%% 6. MATLAB's dtw() for comparison
[dist_matlab, ix, iy] = dtw(x, y);

subplot(3,1,3); hold on; grid on;
plot(tx(ix), x(ix), 'co-', 'MarkerFaceColor','c', 'LineWidth',1.5);
plot(tx(ix), y(iy), 'ms-', 'MarkerFaceColor','m', 'LineWidth',1.8);
xlabel('Time (rad)');
title(sprintf('MATLAB dtw() - Distance = %.3f', dist_matlab));
legend('Warped x', 'Warped y');

%% 7. Cost Matrix + Path
figure(2); clf;
imagesc(ty, tx, D); colormap(flipud(hot)); colorbar;
hold on;
plot(ty(q), tx(p), 'g-', 'LineWidth', 4);          
plot(ty(q), tx(p), 'wo', 'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('y index (sine)'); ylabel('x index (chirp)');
title('DTW Cost Matrix + Optimal Warping Path');
axis xy tight;

%% 8. 3D Visualization: Cost Surface + Optimal Warping Path (BEST ONE)
figure(4); clf;

% Create meshgrid for surf
[Y_grid, X_grid] = meshgrid(ty, tx);   % Note: rows = tx (x), columns = ty (y)

% Plot the accumulated cost as a colored surface
surf(Y_grid, X_grid, D, 'EdgeColor','none');   % or 'EdgeColor','interp' for smoother
shading interp;
colormap(hot(256));        % or parula, jet, etc.
colorbar;
view(90, -90);             % Top-down view (like imagesc but in 3D)
rotate3d on;
hold on;

%% Plot the optimal warping path as a thick green 3D line ON the surface
% We already have p (indices in x) and q (indices in y) from backtracking

% Convert indices → real time values
time_y = ty(q);      % time values for signal y along the path
time_x = tx(p);      % time values for signal x along the path
cost_along_path = D(sub2ind(size(D), p, q));  % actual cost values at each point

% 3D line: (time_y, time_x, cost)
plot3(time_y, time_x, cost_along_path, ...
      'g-', 'LineWidth', 5, 'Color', [0 1 0]);   % Bright green

% Optional: add white markers with black edge for extra clarity
plot3(time_y, time_x, cost_along_path, ...
      'wo', 'MarkerSize', 2, 'MarkerFaceColor','w', 'LineWidth', 1.5);

%% Beautify
xlabel('Signal y time (sine wave) →');
ylabel('Signal x time (chirp) ←');
zlabel('Accumulated DTW Cost');
title('DTW Cost Surface + Optimal Warping Path (Green Line)', 'FontWeight','bold');

% Make axes tight and pretty
axis tight;
grid off;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% Optional: slightly tilt for better 3D feel 
% view(45, 30);
%% Results
fprintf('=== DTW Summary ===\n');
fprintf('Manual total cost      : %.4f\n', D(end,end));
fprintf('MATLAB dtw distance    : %.4f\n', dist_matlab);
fprintf('Warping path length    : %d points\n', length(p));