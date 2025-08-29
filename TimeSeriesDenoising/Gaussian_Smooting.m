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

%% Create a Gaussian kernel

% Full width half-maximum: the key gaussian parameter
fwhm = 25; % = w in ms

% normalized time factor in ms
k = 100;
gtime = 1000*(-k:k)/srate; % = t

% create gaussian window
%           $g=e^{\frac{-4 \ln (2) t^2}{w^2}}$
%           g=e^((-4ln(2)t^(2))/(w^(2)))
gauswin = exp(-(4 * log(2)*gtime.^2 / fwhm^2));


% compute empirical FWHM
PrePeakHalf = k + dsearchn (gauswin(k+1:end)',.5);
PstPeakHalf = dsearchn (gauswin(1:k)',.5);
    empFWHM = gtime(PrePeakHalf)-gtime(PstPeakHalf);

% Show the gaussian 
figure(1), clf, hold on;
plot(gtime, gauswin,"ko-", "MarkerFaceColor","w", "linew", 2);
plot (gtime([PrePeakHalf PstPeakHalf]), gauswin([PrePeakHalf PstPeakHalf]), "m","linew", 1);


% then normalize Gaussian to unit energy
gauswin = gauswin / sum(gauswin);
title([ 'Gaussian kernel with requeted FWHM ' num2str(fwhm) ' ms (' num2str(empFWHM) ' ms achieved)' ])
xlabel('Time (ms)'), ylabel('Gain')

%% implement the filter
%      $y_t=\sum_{i=t-k}^{t+k} x_i g_i$
%      y_(t)=sum_(i=t-k)^(t+k)x_(i)g_(i)

% initialize filtered signal vector
filtsigG = signal;

% implement the running mean filter
for i=k+1:n-k-1
    % each point is the weighted average of k surrounding points
    
    filtsigG(i) = sum( signal(i-k:i+k).*gauswin );
end

% plot
figure(2), clf, hold on
plot(time,signal,'r')
plot(time,filtsigG,'k','linew',3)

xlabel('Time (s)'), ylabel('amp. (a.u.)')
legend({'Original signal';'Gaussian-filtered'})
title('Gaussian smoothing filter')

%% for comparison, plot mean smoothing filter

% initialize filtered signal vector
filtsigMean = zeros(size(signal));

% implement the running mean filter
k = 20; % filter window is actually k*2+1
for i=k+1:n-k-1
    % each point is the average of k surrounding points
    filtsigMean(i) = mean(signal(i-k:i+k));
end

plot(time,filtsigMean,'b','linew',2)
legend({'Original signal';'Gaussian-filtered';'Running mean'})
zoom on

%% done.




