clc;
clear all;
warning off;
close all;

% Input signals 
s=wavread('01 Near Speech');
n=wavread('01 Far Speech(E)');

% Initialization
N=240282;
p=2;
mu=0.02;
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
%figure(1)
 %plot (d)
%sound (d,8000);
%pause(30);

% Algorithm
for i=p:N
    xvec=x(i:-1:i-p+1);
    y(i)=sum(w'*xvec);
    e(i)=d(i)-y(i);
    %den=1+(e(i)*e(i)*beta);
    %mu=beta/den;
    w=w+2*mu*e(i)*xvec;
end
figure(4);
plot(e);
 %sound(e,8000);
% Calculating MSE
for i=1:N
    err(i)=(s(i)-e(i)).^2;
    %nn(i)=n(i).^2;
    %ss(i)=s(i);
end
MSE=mean(err)

% Calculating SNR
% rms_signal= sqrt(mean(e.^2));
% rms_noise= sqrt(mean((ss-e).^2));
% Lsig= 20*log10(rms_signal);
% Lnoise= 20*log10(rms_noise);
% SNR= Lsig - Lnoise