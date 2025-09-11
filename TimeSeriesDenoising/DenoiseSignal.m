
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
supertresh = find(-thershold<origSignal<thershold);
filtsigM = origSignal;
n = length(origSignal);
k= 130;

for ti=1:length(supertresh)
    lowbnd = max(1,supertresh(ti)-k);
    uppbnd = min(supertresh(ti)+k,n);
    filtsigM (supertresh(ti)) = median(origSignal(lowbnd:uppbnd));
end
figure(2),%hold on;
 plot(filtsigM);
 title("Median filter")


%%
% Mean Smooth
k=3
flitsigMean = zeros(size(filtsigM));
for ims=1:length(filtsigM)
    filtsigMean(ims) = mean(filtsigM(max(1,ims-k):min(length(filtsigM),ims+k)));
end
figure(4);
plot(filtsigMean);
title("Meanfilter");


figure(10);
plot(cleanedSignal)