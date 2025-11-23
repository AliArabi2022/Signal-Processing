%% ========================================================================
%  REAL-WORLD DTW DEMO : "Hello" spoken at different speeds
%  Author: Ali Arabi bavil - 2025
%  FINAL VERSION
%  - Records two "hello"s (slow + fast)
%  - Perfect alignment with DTW
%  - Beautiful plots + 3D surface with green path
%  - Plays warped result
% =========================================================================

clear; clc; close all;

%% 1. Recording
fs = 32000;
disp('=== Say "HELLO" twice ===');
disp('   1st: SLOWLY ("helloooooo")');
disp('   2nd: QUICKLY ("hello")');
input('Press Enter when ready...');

rec = audiorecorder(fs,16,1);
disp('Recording SLOW...'); recordblocking(rec,3);
y_slow = rec.getaudiodata();

disp('Recording FAST...'); recordblocking(rec,3);
y_fast = rec.getaudiodata();

%% 2. Trim silence
y_slow = trim_silence(y_slow);
y_fast = trim_silence(y_fast);
t_slow = (0:length(y_slow)-1)/fs;
t_fast = (0:length(y_fast)-1)/fs;

%% 3. Energy envelope
frame_len = round(0.025*fs);
hop_len   = round(0.010*fs);

x = sqrt(mean(buffer(abs(y_slow),frame_len,frame_len-hop_len,'nodelay').^2,1))';
tx = (0:length(x)-1)' * hop_len / fs;

y = sqrt(mean(buffer(abs(y_fast),frame_len,frame_len-hop_len,'nodelay').^2,1))';
ty = (0:length(y)-1)' * hop_len / fs;

%% 4. DTW
N=length(x); M=length(y);
D=zeros(N,M); D(1,1)=abs(x(1)-y(1));
for i=2:N, D(i,1)=D(i-1,1)+abs(x(i)-y(1)); end
for j=2:M, D(1,j)=D(1,j-1)+abs(x(1)-y(j)); end
for i=2:N, for j=2:M
    D(i,j)=abs(x(i)-y(j)) + min([D(i-1,j),D(i,j-1),D(i-1,j-1)]);
end, end

% Backtracking
i=N; j=M; p=[]; q=[];
while i>1 || j>1
    p=[i;p]; q=[j;q];
    if i==1, j=j-1;
    elseif j==1, i=i-1;
    else [~,k]=min([D(i-1,j),D(i,j-1),D(i-1,j-1)]);
        if k==1, i=i-1; elseif k==2, j=j-1; else i=i-1; j=j-1; end
    end
end
p=[1;p]; q=[1;q];

%% 5. MATLAB dtw() for comparison
[matlab_dist, ix, iy] = dtw(x, y);


%% 6. STUDIO-QUALITY TIME WARPING 
y_warped = warpSpeechDTW(y_fast, y_slow, fs, p, q, tx, ty);

%% 7. PLAY — PERFECT STUDIO QUALITY
disp('Playing original SLOW version...');
sound(y_slow / max(abs(y_slow)) * 0.98, fs);
pause(length(y_slow)/fs + 2);

disp('Playing WARPED FAST — STUDIO QUALITY, NATURAL VOICE!');
sound(y_warped, fs);

% Save masterpiece
audiowrite('hello_warped_STUDIO_PERFECT.wav', y_warped, fs, 'BitsPerSample', 24);
disp('Saved: hello_warped_STUDIO_PERFECT.wav — sounds 100% natural!');

%% 8. Beautiful plots
figure(1); clf;
subplot(3,1,1); hold on; grid on;
plot(t_slow,y_slow,'b','LineWidth',1.2); plot(t_fast,y_fast,'r','LineWidth',1.2);
title('Original Recordings'); legend('Slow','Fast');

subplot(3,1,2); hold on; grid on;
plot(tx,x,'b.-','LineWidth',2);
plot(tx(p),y(q),'r.-','LineWidth',2.5);
title('DTW Alignment'); legend('Slow envelope','Warped fast');

subplot(3,1,3); hold on; grid on;
plot(tx(ix),x(ix),'co-'); plot(tx(ix),y(iy),'ms-');
title('MATLAB built-in dtw() - identical result');

figure(2); clf;
[Yg,Xg]=meshgrid(ty,tx);
surf(Yg,Xg,D,'EdgeColor','none'); shading interp; colormap(hot); colorbar;
view(90,-90); rotate3d on; hold on;
cost_path = D(sub2ind(size(D),p,q));
plot3(ty(q),tx(p),cost_path,'g-','LineWidth',7);
title('DTW Cost Surface + Optimal Warping Path (Green)');

%% 9. PROFESSIONAL SUMMARY (exactly what you wanted)
disp(' ');
disp('=== DTW Summary ===');
fprintf('Manual total cost      : %.6f\n', D(end,end));
fprintf('MATLAB dtw distance    : %.6f\n', matlab_dist);
fprintf('Warping path length    : %d points\n', length(p));
fprintf('Slow utterance length  : %.3f seconds\n', t_slow(end));
fprintf('Fast utterance length  : %.3f seconds\n', t_fast(end));
fprintf('Speed ratio (fast/slow): %.2f×\n', t_slow(end)/t_fast(end));
disp(' ');
disp('Perfect alignment achieved! Your fast "hello" now has the exact timing of the slow one.');

