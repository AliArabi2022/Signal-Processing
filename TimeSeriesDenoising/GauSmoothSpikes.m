%% generate time series of random spikes

% number of spikes
n = 300;

% inter-spike intervals (exponential distribution for bursts)
isi = round(exp( randn(n,1) )*10);


figure(1);
subplot(2,1,1);
plot(isi);
subplot(2,1,2);
histogram(isi,100);


% generate spike time series
spikets = 0;
for i=1:n
    spikets(sum(isi(1:i))) = 1;
end


% plot
figure(2), clf, hold on

h = plot(spikets);
set(gca,'ylim',[0 1.01],'xlim',[0 length(spikets)+1])
set(h,'color',[1 1 1]*.7)
xlabel('Time (a.u.)')

%% create and implement Gaussian window

% full-width half-maximum: the key Gaussian parameter
fwhm = 15; % in points

% normalized time vector in ms
k = 100;
gtime = -k:k;

% create Gaussian window
gauswin = exp( -(4*log(2)*gtime.^2) / fwhm^2 );
gauswin = gauswin / sum(gauswin);

% initialize filtered signal vector
filtsigG = zeros(size(spikets));

% implement the weighted running mean filter
for i=k+1:length(spikets)-k-1
    filtsigG(i) = sum( spikets(i-k:i+k).*gauswin );
end


% plot the filtered signal (spike probability density)
plot(filtsigG,'r','linew',2)
legend({'Spikes','Spike p.d.'})
title('Spikes and spike probability density')
zoom on

%% done.
