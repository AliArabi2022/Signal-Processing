clc
clear all
%% Start
srate = 1000;        % نرخ نمونه‌برداری
time  = -4:1/srate:4; % بازه‌ی زمانی از -4 تا 4 ثانیه
pnts  = length(time); % تعداد نقاط
% ساخت فیلتر سینک
f = 8;
sincfilt = sin(2*pi*f*time) ./ time;
% adjust NaN and normalize filter to unit-gain
sincfilt(~isfinite(sincfilt)) = max(sincfilt);
% چون در
% t=0 
% مقدارش 
% NaN
sincfilt = sincfilt./sum(sincfilt);
%نرمالایزش می‌کنیم
sincfiltW = sincfilt .* hann(pnts)';
% جای اینکه سینک بی‌نهایت باشه (از -∞ تا +∞)
%  اون رو برش می‌زنیم و با یک پنجره هموار مثل
% Hann 
% ضرب می‌کنیم
%  این باعث میشه در حوزه‌ی فرکانس لوب‌های جانبی
% (side lobes)
% کاهش پیدا کنن → پاسخ فرکانسی تمیزتر

%%% plotting
figure(1), clf

% sinc filter
subplot(221)
plot(time,sincfilt,'k','linew',1)
xlabel('Time (s)')
title('Non-windowed sinc function')

% now plot the windowed sinc filter
subplot(222)
plot(time,sincfiltW,'k','linew',1)
xlabel('Time (s)')
title('Windowed sinc function')

% sinc filter
subplot(223), hold on
hz = linspace(0,srate/2,floor(pnts/2)+1);
pw = abs(fft(sincfilt));
plot(hz,pw(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 f*3],'YScale','lo','ylim',[10e-7 10])
plot([1 1]*f,get(gca,'ylim'),'r--')
xlabel('Frequency (Hz)'), ylabel('Gain')
% now plot the windowed sinc filter
subplot(224), hold on
hz = linspace(0,srate/2,floor(pnts/2)+1);
pw = abs(fft(sincfiltW));
plot(hz,pw(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 f*3],'YScale','lo','ylim',[10e-7 10])
plot([1 1]*f,get(gca,'ylim'),'r--')
xlabel('Frequency (Hz)'), ylabel('Gain')



%% apply the filter to noise

% generate data as integrated noise
data = cumsum( randn(pnts,1) );

% reflection
datacat = [data; data(end:-1:1)];

% apply filter (zero-phase-shift)
dataf = filter(sincfiltW,1,datacat);
dataf = filter(sincfiltW,1,dataf(end:-1:1));

% flip forwards and remove reflected points
dataf = dataf(end:-1:pnts+1);


% compute spectra of original and filtered signals
hz = linspace(0,srate/2,floor(pnts/2)+1);
powOrig = abs(fft(data)/pnts).^2;
powFilt = abs(fft(dataf)/pnts).^2;


% plot
figure(2), clf
subplot(211)
plot(time,data, time,dataf,'linew',1)
xlabel('Time (s)'), ylabel('Amplitude')
legend({'Original signal';'Windowed-sinc filtered'})
zoom on


% plot original and filtered spectra
subplot(212)
plot(hz,powOrig(1:length(hz)), hz,powFilt(1:length(hz)), 'linew',1)
set(gca,'xlim',[0 40],'YScale','log')
xlabel('Frequency (Hz)'), ylabel('Power')
legend({'Original signal';'Windowed-sinc filtered'})
%% with different windowing functions

sincfiltW = cell(3,1);

tapernames = {'Hann';'Hamming';'Gauss'};

% with Hann taper
%sincfiltW{1} = sincfilt .* hann(pnts)';
hannw = .5 - cos(2*pi*linspace(0,1,pnts))./2;
sincfiltW{1} = sincfilt .* hannw;


% with Hamming taper
%sincfiltW{2} = sincfilt .* hamming(pnts)';
hammingw = .54 - .46*cos(2*pi*linspace(0,1,pnts));
sincfiltW{2} = sincfilt .* hammingw;


% with Gaussian taper
sincfiltW{3} = sincfilt .* exp(-time.^2);



figure(3), clf

for filti=1:length(sincfiltW)
    subplot(211), hold on
    plot(time,sincfiltW{filti},'linew',2)
    xlabel('Time (s)'), title('Time domain')
    
    subplot(212), hold on
    hz = linspace(0,srate/2,floor(pnts/2)+1);
    pw = abs(fft(sincfiltW{filti}));
    plot(hz,pw(1:length(hz)),'linew',3)
    set(gca,'xlim',[f-3 f+10],'YScale','lo','ylim',[10e-7 10])
    xlabel('Frequency (Hz)'), ylabel('Gain')
    title('Frequency domain')
    
end

legend(tapernames)

%% done.

