%% ایجاد سیگنال اولیه
% یک بخش صفر + یک موج سینوسی کوتاه + دوباره بخش صفر
signal = [ zeros(1,200) sin(linspace(0,4*pi/2,50)) zeros(1,200) ];
n = length(signal);   % طول سیگنال

%% رسم سیگنال و طیف فرکانسی آن
figure(1), clf
subplot(221)
plot(1:n,signal,'ko-')          % رسم سیگنال در حوزه زمان
set(gca,'xlim',[0 n+1])
title('Original signal')
xlabel('Time points (a.u.)')

subplot(222)
plot(linspace(0,1,n),abs(fft(signal)),'ko-','markerfacecolor','w')  % طیف فرکانس
set(gca,'xlim',[0 .5])
xlabel('Frequency (norm.)'), ylabel('Energy')
title('Frequency-domain signal representation')

%% طراحی و اعمال فیلتر پایین‌گذر
% fir1: طراحی فیلتر FIR با طول 50 و فرکانس قطع 0.6 (به نسبت نایکوئیست)
fkern = fir1(50,.6,'low');  

% فیلتر کردن سیگنال به صورت علی (forward)
fsignal = filter(fkern,1,signal);
subplot(234)
plot(1:n,signal, 1:n,fsignal)
set(gca,'xlim',[0 n+1]), axis square
xlabel('Time (a.u.)')
legend({'Original';'Forward filtered'})

% نمایش طیف سیگنال فیلترشده روی نمودار فرکانس
fsignalFlip = fsignal(end:-1:1);   % برگرداندن سیگنال به عقب
subplot(222), hold on
plot(linspace(0,1,n),abs(fft(fsignal)),'r','linew',3)

%% فیلتر کردن سیگنال برعکس (Backward filtering)
fsignalF = filter(fkern,1,fsignalFlip);
subplot(235)
plot(1:n,signal, 1:n,fsignalF)
set(gca,'xlim',[0 n+1]), axis square
xlabel('Time (a.u.)')
legend({'Original';'Backward filtered'})

%% بازگرداندن سیگنال به ترتیب صحیح و رسم خروجی نهایی (Zero-phase)
fsignalF = fsignalF(end:-1:1);   % دوباره برعکس کردن (جلو + عقب فیلتر شده)
subplot(236)
plot(1:n,signal, 1:n,fsignalF)
set(gca,'xlim',[0 n+1]), axis square
xlabel('Time (a.u.)')
legend({'Original';'Zero-phase filtered'})

% نمایش طیف سیگنال بدون شیفت فاز
subplot(222)
plot(linspace(0,1,n),abs(fft(fsignalF)),'mo','markersize',7,'markerfacecolor','w','linew',2)
legend({'Original';'Forward filtered';'Zero-phase-shift'})

%% پایان
