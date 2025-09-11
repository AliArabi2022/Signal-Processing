clc;
clear all;
%%
load Data\denoising_codeChallenge.mat

%%
% Meadian Filter
% figure(1), clf
% histogram(origSignal,100);
% zoom on

thershold = 4;
%positive 
supertresh = find(origSignal>thershold);
filtsigM = origSignal;
n = length(origSignal);
k= 5;

for ti=1:length(supertresh)
    lowbnd = max(1,supertresh(ti)-k);
    uppbnd = min(supertresh(ti)+k,n);
    filtsigM (supertresh(ti)) = median(origSignal(lowbnd:uppbnd));
end


supertresh0 = find(origSignal<-thershold);
for ti=1:length(supertresh0)
    lowbnd = max(1,supertresh0(ti)-k);
    uppbnd = min(supertresh0(ti)+k,n);
    filtsigM (supertresh0(ti)) = median(origSignal(lowbnd:uppbnd));
end


figure(2),%hold on;
plot(filtsigM);
title("Median filter")

k2=150
flitsigMean = zeros(size(filtsigM));
for ims=1:length(filtsigM)
    filtsigMean(ims) = mean(filtsigM(max(1,ims-k2):min(length(filtsigM),ims+k2)));
end
figure(4);subplot(211)
plot(1:length(filtsigM),filtsigMean);
title("Meanfilter");
subplot(212)
plot(1:length(filtsigM),cleanedSignal);
title("cleanedSignal");
