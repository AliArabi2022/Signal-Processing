%%
clc
clear
% close all



%%

% Create signal
srate = 1000; % Hz sample rate
time = 0:1/srate:3; % splite from 1 to 3| step = 1/srate
n = length(time); % number of sample
p = 15; % poles for random interpolation

% Noise level, measured in standard deviation
noiseamp = 5;

% amplitude modulator and noise level
ampl = interp1(rand(p,1)*30, linspace(1,p,n));
noise = noiseamp * randn(size(time));
signal = ampl + noise;
 subplot(2,1,1);
 plot(time,ampl);
 subplot(2,1,2);
 plot(time,signal);
 pause; % Pauses until a key is pressed
 disp('Key pressed, continuing.');

% initialize filtered signal vector
filtsig = zeros(size(signal));

% implement the running mean filter
k = 20; % filter window is k*2+1

for i=k+1:n-k-1
    % each point is the avarage of k surrounding points
    filtsig(i)= mean(signal(i-k:i+k));
end

% compute the window size in ms
windowsize = 1000*(k*2+1)/ srate;

% plot the noisy and filtered signals
figure(1), clf, hold on
plot(time, signal, time, filtsig, "linew", 2)

% draw a patch to indicate the window size
tidx = dsearchn(time', 1);
ylim = get (gca, "ylim");
patch(time([tidx-k tidx-k tidx-k tidx-k ]),ylim([1 2 2 1]), "k", ...
    "facealpha", .25, "linestyle", "none")
plot (time([tidx tidx]), ylim, "k--")

xlabel("Time(sec.)"), ylabel("Amplitude")
title([ "Running-mean filter with a k=" num2str(round(windowsize)) "-ms filter"])
legend({"Signal" ; "Filtered" ; "Window" ; "Window Center"})
zoom on

%% done

