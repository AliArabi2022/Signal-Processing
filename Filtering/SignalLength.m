% parameters
dataN = 10000;   % طول سیگنال
filtN = 5001;    % مرتبه فیلتر FIR (خیلی بزرگ!)

% generate data
signal = randn(dataN,1);  % سیگنال نویزی تصادفی

% create filter kernel
fkern = fir1(filtN,.01,'low');  % فیلتر FIR پایین‌گذر با cutoff=0.01 نایکوئیست

% apply filter kernel to data
% fdata = filtfilt(fkern,1,signal);


% use reflection to increase signal length!
signalRefl = [ signal(end:-1:1); signal; signal(end:-1:1) ];

% apply filter kernel to data
fdataR = filtfilt(fkern,1,signalRefl);

% and cut off edges
fdata = fdataR(dataN+1:end-dataN);


%%
% parameters
dataN = 10000;   % طول سیگنال
filtN = 5001;    % مرتبه فیلتر FIR (خیلی بزرگ!)

% generate data
signal = randn(dataN,1);  % سیگنال نویزی تصادفی

% create filter kernel
fkern = fir1(filtN,.01,'low');  % فیلتر FIR پایین‌گذر با cutoff=0.01 نایکوئیست


% fkern فرض بر این است که از قبل ساخته شده (مثلاً fkern = fir1(filtN,...))

lenb = length(fkern);              % تعداد ضرایب b
minLen = 3 * lenb + 1;            % حداقل طول لازم (strict > در عمل +1 برای عدد صحیح)
dataN = length(signal);

if dataN < minLen
    pad = minLen - dataN;         % تعداد نمونه‌ای که باید اضافه کنیم
    % نمونه‌های لازم از ابتدا و انتها را به صورت بازتاب برداریم
    left  = signal(pad:-1:1);                % pad نمونه از ابتدا، معکوس
    right = signal(end:-1:end-pad+1);        % pad نمونه از انتها، معکوس
    signalRefl = [left; signal; right];
else
    pad = 0;
    signalRefl = signal;
end

% حالا 
% filtfilt
% را روی سیگنال بازتاب‌شده اجرا کن
fdataRefl = filtfilt(fkern, 1, signalRefl);

% برش دادن بخش اصلی
fdata = fdataRefl(pad+1 : pad+dataN);

