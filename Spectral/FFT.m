
clc;
clear all;

%% Generate a multispectral noisy signal

% simulation parameters
Srate = 1234; % in Hz
npnts = Srate*2;
time = (0:npnts-1)/Srate;
% frequencies to include
frex = [ 12, 18, 30];
signal = zeros (size(time));
% loop over frequencies to create signal
for fi = 1:length(frex)
    signal = signal + fi*sin(2*pi*frex(fi)*time);
end

signal = signal +5 * sin(2*pi*36*time);
% add some noise 
signal = signal + randn(size(signal));

%amplitude spectrum via Fourier Transform
signalX = fft(signal);
signalAmp = 2*abs(signalX)/npnts;
% vector of frequencies in Hz
hz = linspace(0,Srate/2, floor(npnts/2)+1);

%% plots

figure(1), clf
subplot(211);
plot(time, signal)
zoom on;
xlabel('Time (s)'), ylabel('Amplitude')
title('Time domain')
subplot(212)
stem(hz,signalAmp(1:length(hz)), "ks", "linew", 2, "markersize",10);
set(gca, "xlim", [0 max(frex)*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')
subplot(211), hold on;
plot(time, ifft(signalX), "ro")
legend({'Original';'IFFT reconstructed'})
%%
% Example with Iran Real data
search = load ("Data\SignalProcessingIran.mat");
searchdata = (search.signalProcessingIran)';
N = length(searchdata);

searchpow = abs(fft(searchdata));
hz = linspace(0,52,N);

figure(2), clf
subplot(211)
plot(searchdata,'ko-','markerfacecolor','m','linew',1)
xlabel('Time (weeks)'), ylabel('Search volume')

subplot(212)
plot(hz,searchpow,'ms-','markerfacecolor','k','linew',1)
xlabel('Frequency (norm.)'), ylabel('Search power')
set(gca,'xlim',[0 12])
