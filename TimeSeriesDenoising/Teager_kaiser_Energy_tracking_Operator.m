%%
% Denoise EMG signals via Teager Kaiser Energy traking Oprator
%
%
% In discrete-time applications, a common form of the operator is:
%   Ψ[s(n)] = s(n)^2 - s(n-1)s(n+1)
%  where s(n) is the signal at the current time point,
% and s(n-1) and s(n+1) are the previous and next time points. 

% import data
load Data\emg4TKEO.mat

% initialize filtered signal
emgf = emg;
% figure(1);
% plot(emgtime,emg);


% the loop version for interpretability
for i=2:length(emgf)-1
    emgf(i) = emg(i)^2 - emg(i-1)*emg(i+1);
end

% the vectorized version for speed and elegance
emgf2 = emg;
emgf2(2:end-1) = emg(2:end-1).^2 - emg(1:end-2).*emg(3:end);

%% convert both signals to zscore

% find timepoint zero
time0 = dsearchn(emgtime',0);

% convert original EMG to z-score from time-zero
emgZ = (emg-mean(emg(1:time0))) / std(emg(1:time0));

% same for filtered EMG energy
emgZf = (emgf-mean(emgf(1:time0))) / std(emgf(1:time0));

%% plot

figure(2), clf

% plot "raw" (normalized to max-1)
subplot(211), hold on
plot(emgtime,emg./max(emg),'b','linew',2)
plot(emgtime,emgf./max(emgf),'m','linew',2)

xlabel('Time (ms)'), ylabel('Amplitude or energy')
legend({'EMG';'EMG energy (TKEO)'})


% plot zscored
subplot(212), hold on
plot(emgtime,emgZ,'b','linew',2)
plot(emgtime,emgZf,'m','linew',2)

xlabel('Time (ms)'), ylabel('Zscore relative to pre-stimulus')
legend({'EMG';'EMG energy (TKEO)'})

%% done.