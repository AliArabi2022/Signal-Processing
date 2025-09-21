clc;
clear all;

%% Notch filter design
fs      = 350; % نرخ نمونه‌برداری (Hz)
timevec = (0:fs*7-1)/fs;  % بردار زمان برای 7 ثانیه
npnts   = length(timevec);

% سیگنال اصلی (random walk = انتگرال نویز سفید)
yOrig = cumsum(randn(npnts,1));

% سیگنال اندازه‌گیری شده: سیگنال اصلی + نویز سفید قوی + سینوسی 50Hz
y = yOrig + 50*randn(npnts,1) + 40*sin(2*pi*50*timevec)';



notchfreq = 50;     % فرکانس مورد نظر برای حذف
bw        = 4;      % پهنای باند notch (Hz)
order     = round(7*fs/notchfreq);

% تعریف شکل فیلتر: 1 همه‌جا، 0 در اطراف 50Hz
shape = [1 1 0 0 1 1];
frex  = [0 notchfreq-bw notchfreq notchfreq+1 notchfreq+bw fs/2]/(fs/2);

% طراحی فیلتر notch
filtkern = firls(order, frex, shape);

filtkernX = abs(fft(filtkern,npnts)).^2;
% اطمینان از اینکه ضرایب ستونی باشن
%filtkern = filtkern(:);

% اعمال فیلتر
yNotch = filtfilt(filtkern,1,y);

%% Plot results
figure, clf
subplot(211)
plot(timevec,y, 'color',[.6 .6 .6])
hold on
plot(timevec,yNotch,'b','linew',1.5)
xlabel('Time (sec)'), ylabel('Amplitude')
legend({'Original noisy signal','After 50Hz Notch'})
title('Time domain')

% Power spectrum
yNotchX = abs(fft(yNotch)/npnts).^2;
hz = linspace(0,fs/2,floor(npnts/2)+1);

subplot(212)
plot(hz,y(1:length(hz)),'k','linew',1), hold on
plot(hz,yNotchX(1:length(hz)),'b','linew',2)
set(gca,'yscale','log','xlim',[0 80])
legend({'Original','Notch-filtered'})
xlabel('Frequency (Hz)'), ylabel('Power (log)')
title('Frequency domain')




figure(2), clf
subplot(321)
plot((-order/2:order/2+1)/fs,filtkern,'k','linew',3)
xlabel('Time (s)')
title('Filter kernel')

subplot(322), hold on
plot(frex*fs/2,shape,'r','linew',1)
plot(hz,filtkernX(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 60])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Filter kernel spectrum')


%%% now apply the filter to the data
yFilt = filtfilt(filtkern,1,y);

subplot(312)
h = plot(timevec,y, timevec,yFilt,'linew',2);
set(h(1),'color',[1 1 1]*.4)
legend({'Signal';'Filtered'})
xlabel('Time (sec.)'), ylabel('Amplitude')


%%% power spectra of original and filtered signal
yOrigX = abs(fft(y)/npnts).^2;
yFiltX = abs(fft(yFilt)/npnts).^2;

subplot(313)
plot(hz,yOrigX(1:length(hz)), hz,yFiltX(1:length(hz)),'linew',2);
set(gca,'xlim',[0 fs/5],'yscale','log')
legend({'Signal';'Filtered'})
xlabel('Frequency (Hz)'), ylabel('Power (log)')


