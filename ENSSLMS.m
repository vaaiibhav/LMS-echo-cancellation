clc;
clear all;
warning off;
close all;

% Input signals 
s=wavread('E:\mru\ME\Seminar n project\database\clean\sp11.wav');
n=wavread('E:\mru\ME\Seminar n project\database\noise signal\airport\0db.wav');

% Initialization
N=22529;
p=2;
mu=2;
w=zeros(p,1);
x=zeros(N,1);
d=zeros(N,1);

% Noise and noisy inputs
for i=1:N
    v1(i)=0.9*n(i);
    v2(i)=(-0.9)*n(i);
    x(i)=v2(i);
    d(i)=s(i)+v1(i);
end

% Algorithm
for i=p:N
    xvec=x(i:-1:i-p+1);
    y(i)=w'*xvec;
    e(i)=d(i)-y(i);
%     den=1+(e(i)*e(i)*beta);
%     mu=beta/den;
    w=w+mu*e(i)*xvec;
end

% Calculating MSE
for i=1:N
    err(i)=(s(i)-e(i)).^2;
    nn(i)=n(i).^2;
    ss(i)=s(i);
end
MSE=mean(err)

% Calculating SNR
rms_signal= sqrt(mean(e.^2));
rms_noise= sqrt(mean((ss-e).^2));
Lsig= 20*log10(rms_signal);
Lnoise= 20*log10(rms_noise);
SNR= Lsig - Lnoise

% Calculating Misadjustment
mMSE=mean(nn);
Mi=(MSE-mMSE)/mMSE
 
% % Output file wavwrite
% wavwrite(e,'E:\mru\ME\Seminar n project\database\Estimated_op\ENLMS\street\15db\sp11');

% Plots for ip/op
% for i =1:1:N
%      S(i)= s(i)';
%      E(i)=e(i);
% end
n= 1:1:N;
plot(n,e,'g');
title('Actual speech signal versus Estimated signal');
xlabel('time index'); ylabel('signal value');
hold on;
plot(n,ss,'r');
% n= 1:1:N;
% plot(n,d,'g');
% title('Comparision of Estimated speech signal with noisy and clean speech signals');
% xlabel('time index'); ylabel('signal value');
% hold on;
% plot(n,e);
% hold on;
% plot(n,s,'r');

% Spectrogram representation
% figure;
% subplot 411;
% [S,F,T,P]=spectrogram(s,512,256,512,44100);
% surf(T,F,10*log10(abs(P)),'edgecolor','none')
% colormap(jet);
% axis tight;
% view(0,90);
% title('Spectrogram of Clean speech signal');
% xlabel('Time');
% ylabel('Hz');
% subplot 412;
% [S,F,T,P]=spectrogram(x,512,256,512,44100);
% surf(T,F,10*log10(abs(P)),'edgecolor','none')
% colormap(jet);
% axis tight;
% view(0,90);
% title('Spectrogram of Noise signal');
% xlabel('Time');
% ylabel('Hz');
% subplot 413;
% [S,F,T,P]=spectrogram(d,512,256,512,44100);
% surf(T,F,10*log10(abs(P)),'edgecolor','none')
% colormap(jet);
% axis tight;
% view(0,90);
% title('Spectrogram of Noisy speech signal');
% xlabel('Time');
% ylabel('Hz');
% subplot 414;
% [S,F,T,P]=spectrogram(e,512,256,512,44100);
% surf(T,F,10*log10(abs(P)),'edgecolor','none')
% colormap(jet);
% axis tight;
% view(0,90);
% title('Spectrogram of Estimated clean speech signal');
% xlabel('Time');
% ylabel('Hz');