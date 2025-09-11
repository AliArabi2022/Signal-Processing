% load dataset
load Data\templateProjection.mat


% initialize residual data
residual = zeros(size(EEGdat));

% loop over trials
for i=1 :size(residual,2)

    % build the least-squares model as intercept and EOG from this trial
    X = [ones(size(residual,1),1),eyedat(:,i)];
    % compute regression coefficients for EEG channel
    beta = (X'*X)\(X'*EEGdat(:,i));
    % predicted data
    yHat = X * beta;
    % new data are the residuals after projecting out the best EKG fit
    residual(:,i) = (EEGdat(:,i) - yHat)';
end
%% plotting

% trial averages
figure(1), clf
plot(timevec,mean(eyedat,2), timevec,mean(EEGdat,2), timevec,mean(residual,2),'linew',2)
legend({'EOG';'EEG';'Residual'})
xlabel('Time (ms)')


% show all trials in a map
clim = [-1 1]*20;

figure(2), clf
subplot(131)
imagesc(timevec,[],eyedat')
set(gca,'clim',clim)
xlabel('Time (ms)'), ylabel('Trials')
title('EOG')


subplot(132)
imagesc(timevec,[],EEGdat')
set(gca,'clim',clim)
xlabel('Time (ms)'), ylabel('Trials')
title('EEG')


subplot(133)
imagesc(timevec,[],residual')
set(gca,'clim',clim)
xlabel('Time (ms)'), ylabel('Trials')
title('Residual')

%% done.
