%% create multichannel signal with multiple sampling rates

% initialize signals, time vectors, and sampling rates
% note: hard-coded to three signals
clc;
clear all;
close all;



[fs,timez,signals] = deal( cell(5,1) );


% sampling rates in Hz
fs{1} = 10;
fs{2} = 40;
fs{3} = 83;
fs{4} = 95;
fs{5} = 108;
% create signals
for si=1:length(fs)
    
    % create signal
    signals{si} = cumsum( sign(randn(fs{si},1)) );
    
    % create time vector
    timez{si} = (0:fs{si}-1)/fs{si};
end



% plot all signals
figure(1), clf, hold on

color = 'kbrgc';
shape = 'os^*+';
for si=1:length(fs)
    plot(timez{si},signals{si},[ color(si) shape(si) '-' ],'linew',1,'markerfacecolor','w','markersize',6)
end
axlims = axis;
xlabel('Time (s)')
legend( 'sig1','sig2','sig3','sig4','sig5')

%% upsample to fastest frequency

% in Hz
[newSrate,whichIsFastest] = max(cell2mat(fs));

% need to round in case it's not exact
newNpnts = round( length(signals{whichIsFastest}) * (newSrate/fs{whichIsFastest}) );

% new time vector after upsampling
newTime = (0:newNpnts-1) / newSrate;

%% continue on to interpolation

% initialize (as matrix!)
newsignals = zeros(length(fs),length(newTime));

for si=1:length(fs)
    
    % define interpolation object
    F = griddedInterpolant(timez{si},signals{si},'spline');
    
    % query that object at requested time points
    newsignals(si,:) = F(newTime);
end

%% plot for comparison


% plot all signals
figure(2), clf, hold on

for si=1:length(fs)
    plot(newTime,newsignals(si,:),[ color(si) shape(si) '-' ],'linew',1,'markerfacecolor','w','markersize',6)
end

% set axis limits to match figure 1
axis(axlims)

for si=1:length(fs)
    plot(newTime,newsignals(si,:),[ color(si) shape(si) '-' ],'linew',1,'markerfacecolor','w','markersize',6)
end
legend( 'upsig1','upsig2','upsig3','upsig4','sig5')
% set axis limits to match figure 1
axis(axlims)
