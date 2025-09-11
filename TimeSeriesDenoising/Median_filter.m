%%
% Median filter to remove spike noise
% Median filter is non-linear filter
% Apply to selected points, not all data points
% Example:
% [0,1,1,2,0,3,900]
% step1: Sort=> [0,0,1,1,2,3,900]
% step2: Median = middle value
% Or average middle values for example [0,1,1,2,7,0,3,900]
% [0,0,1,1,2,3,7,900]
% (1+2)/2=1.5 
%
%%
clc;
clear all;


% create signal
n = 2000;
signal = cumsum(randn(n,1));


% proportion of time points to replace with noise
propnoise = .05;

% find noise points
noisepnts = randperm(n); 
% randperm(n)
% returns a row vector containing
% a random permutation of the integers from 1 to n without
% repeating elements



noisepnts = noisepnts(1:round(n*propnoise));
% pick the noisepnts's elements from 1 to round(n*propnoise)




% generate signal and replace points with noise
signal(noisepnts) = 50 + rand(size(noisepnts))*100;
% replace the value of signal(noisepnts)
% with 
% amplitudes = rand * 100 + 50


% use hist to pick threshold
figure(1), clf
hist(signal,100);
zoom on

prompt = "What is the threshold value? ";
th = input(prompt)
% visual-picked threshold
threshold = th;


% find data values above the threshold
suprathresh = find( signal>threshold );
% suprathresh = 1:n, % If use median filter for all n point of signal

% initialize filtered signal
filtsig = signal;

% loop through suprathreshold points and set to median of k
k = 20; % actual window is k*2+1
for ti=1:length(suprathresh)
    
    % find lower and upper bounds
    lowbnd = max(1,suprathresh(ti)-k);
    uppbnd = min(suprathresh(ti)+k,n);
    
    % compute median of surrounding points
    filtsig(suprathresh(ti)) = median(signal(lowbnd:uppbnd));
end

% plot
figure(2), clf
plot(1:n,signal, 1:n,filtsig, 'linew',2)
zoom on

%% done.

