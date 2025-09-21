%% پاک کردن محیط کاری
clear

%% ایجاد یک سیگنال تصادفی باند-پَس با نویز
N  = 500;                       % طول سیگنال (تعداد نمونه‌ها)
hz = linspace(0,1,N);           % محور فرکانس نرمال‌شده (0 تا 1= نایکوئیست)

% ساخت یک
% envelope
% گوسی در حوزه فرکانس، با مرکز حدود 0.1 نایکوئیست
gx = exp( -(4*log(2)*(hz-.1)/.1).^2 )*N/2;

% تولید سیگنال در حوزه زمان:
% 1- ضرب 
% envelope
% گوسی در فازهای تصادفی در فرکانس
% 2- ifft → برگشت به حوزه زمان
% 3- اضافه کردن نویز گاوسی (randn)
signal = real(ifft( gx.*exp(1i*rand(1,N)*2*pi) )) + randn(1,N);

%% نمایش سیگنال و طیف توان آن
figure(1), clf

subplot(311)
plot(1:N,signal,'k')
set(gca,'xlim',[0 N+1])
title('Original signal')
xlabel('Time points (a.u.)')

subplot(324)
plot(hz,abs(fft(signal)).^2,'k')
set(gca,'xlim',[0 .5])
xlabel('Frequency (norm.)'), ylabel('Energy')
title('Frequency-domain signal representation')

%% طراحی و اعمال فیلتر پایین‌گذر بدون reflection
order = 150;                            % مرتبه فیلتر FIR
fkern = fir1(order,.6,'low');           % فیلتر پایین‌گذر با fir1

% اعمال فیلتر دوبل (forward-backward) برای حذف شیفت فاز
fsignal = filter(fkern,1,signal);          % فیلتر در جهت مستقیم
fsignal = filter(fkern,1,fsignal(end:-1:1)); % فیلتر در جهت معکوس
fsignal = fsignal(end:-1:1);               % برگرداندن به ترتیب اصلی

% نمایش سیگنال اصلی و فیلتر شده (بدون reflection)
subplot(323), hold on
plot(1:N,signal,'k')
plot(1:N,fsignal,'m')
set(gca,'xlim',[0 N+1])
xlabel('Time (a.u.)')
title('Time domain')
legend({'Original';'Filtered, no reflection'})

% نمایش طیف توان سیگنال فیلتر شده
subplot(324), hold on
plot(hz,abs(fft(fsignal)).^2,'m')
title('Frequency domain')
legend({'Original';'Filtered, no reflection'})

%% اعمال فیلتر با reflection برای کاهش اثر لبه‌ها
% بازتاب دادن ابتدا و انتهای سیگنال به اندازه مرتبه فیلتر
reflectsig = [ signal(order:-1:1) signal signal(end:-1:end-order+1) ];

% اعمال فیلتر دوبل روی سیگنال بازتاب‌دار
reflectsig = filter(fkern,1,reflectsig);           
reflectsig = filter(fkern,1,reflectsig(end:-1:1)); 
reflectsig = reflectsig(end:-1:1);                

% حذف بخش بازتابی و بازگرداندن سیگنال به طول اصلی
fsignal = reflectsig(order+1:end-order);

% نمایش سیگنال اصلی و فیلتر شده (با reflection)
subplot(325), hold on
plot(1:N,signal,'k')
plot(1:N,fsignal,'m')
set(gca,'xlim',[0 N+1])
xlabel('Time (a.u.)')
title('Time domain')
legend({'Original';'Filtered, with reflection'})

% نمایش طیف توان مقایسه‌ای
subplot(326), hold on
plot(hz,abs(fft(signal)).^2,'k')
plot(hz,abs(fft(fsignal)).^2,'m')
legend({'Original';'Filtered, with reflection'})
set(gca,'xlim',[0 .5])
xlabel('Frequency (norm.)'), ylabel('Energy')
title('Frequency domain')

%% نکته: دستور زیر (filtfilt) کار همه این مراحل را به صورت آماده انجام می‌دهد
% fsignal1 = filtfilt(fkern,1,signal);
