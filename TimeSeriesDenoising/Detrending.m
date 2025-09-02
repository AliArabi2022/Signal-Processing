%%





%%
% Create signal with linear trend imposed

n = 2000; 
signal = cumsum(randn(1,n))+ linspace(-30,30,n);

% Linear detrending
detsignal = detrend(signal);

% Plot signal and detrended signal
figure(), clf
plot(1:n,signal,1:n, detsignal,"linew",.3)
legend({ ['Orginal (mean=' num2str(mean(signal)) ')'];['Detrend (mean=' num2str(mean(detsignal)) ')'] });