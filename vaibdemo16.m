clc;
clear all;
close all;
%--------------------------------------------------------------------- 
%Load Data 
[Fname1,Pname1] = uigetfile('01 Near Speech.wav','Select nearspeech FIle'); 
[v,Fs1] = audioread(strcat(Pname1,Fname1));  %Near-end signal  


[Fname2,Pname2] = uigetfile('01 Far speech.wav','Select farspeech FIle'); 
[x,Fs2] = audioread(strcat(Pname2,Fname2));    %Far-end signal 


[Fname3,Pname3] = uigetfile('RIR.wav','Select RIR FIle'); 
[h,Fs3] = audioread(strcat(Pname3,Fname3));    %Far-end signal 

% [h, Fs, nbits] = wavread('c:/audiofiles/room_impulse_response_128taps'); %Room impulse response 
  
%Declare the needed variables  
L=128;        %Length of adaptive filter (same length of RIR)
N=length(x);        %Number of iterations 
lambda_DTD=0.95;    %Constant for calculating decision statistic of DTD 
% DTDbegin=21000;     %The time to activate DTD 
  
%Intial value 0 
w=zeros(L,1);       %Initial weight vector of AF Lx1 
xin=zeros(L,1);     %Initial input signal of AF Lx1 
varMIC=zeros(N,1);  %Initial variance of microphone signal of AF Nx1 
r_em=zeros(N,1);    %Initial Cross correlation between error and microphone signals 
  
%Ambient noise 
WhiteNoise = wgn(N,1,-65);     %With make SNR of 45dB 
%Microphone signal 
EchoSignal=filter(h,1,x);      %Echo signal after filter H 
d=EchoSignal+WhiteNoise+v;     %Desired signal (Microphone Signal) 
 sound(d,Fs1); 
%Make column vectors 
x=x(:);             %Far end signal Nx1 
d=d(:);             %Desired signal Nx1 
  
%The values for calculate Step-Size of Adaptive Filter           
mu=0.014; 
            
%Calculate the average SNR (desired signal/noise) 
powerMic = sum(abs(d).^2)/N;            %Power of Microphone signal 
powerN = sum(abs(WhiteNoise).^2)/N;     %Power of White Noise 
SNR=10*log10(powerMic/powerN);          %Calculate the SNR 
  
%---------------------------------------------------------------------- 
%-------------LMS algorithm for Adaptive Filter----------------------- 
for i=1:N 

    xin(1)=x(i);               %Insert new sample at beginning of input 
     
    y(i)=w'*xin;               %Output signal after adaptive filter    
    error=d(i)-y(i);           %ERROR  
    e(i)=error;                %Store estimation error 
    wtemp = w + 2*mu*error*xin;%Update filter   
 
%-------------ERLE------------------------------------- 
powerD(i) = abs(d(i))^2;    %Power of Microphone signal 
powerE(i)=abs(e(i))^2;      %power of Error signal 
%--------------MSE------------------------------------- 
mse_iteration(i)=error^2;  %Square Error 
end 
for i=1:N-L 
    %MSE - Mean Square Error 
    mse(i)=mean(mse_iteration(i:i+L)); 
    %Echo Return Loss Enhancement 
    ERLE(i)=10*log10(mean(powerD(i:i+L))/mean(powerE(i:i+L))); 
     
     
end 
%---------------------------------------------------------------------- 
%PlOTTING THE NECESSARY SIGNALS 
%---------------------------------------------------------------------- 
figure(1) 
%-------echo signal------------------------- 
subplot(4,1,1) 
plot(EchoSignal) 
xlabel('time (samples)');  
ylabel('echo(n)'); 
title('ECHO SIGNAL: echo(n)') 
grid on 
axis([0 N -1 1]); 
  
%-------Desired signal----------------------- 
subplot(4,1,2) 
plot(d) 
xlabel('time (samples)');  
ylabel('d(n)'); 
title('DESIRED SIGNAL: d(n)') 
grid on 
axis([0 N -1 1]); 
%-------Output signal x(n)------------------- 
subplot(4,1,3) 
plot(y) 
xlabel('time (samples)');  
ylabel('y(n)'); 
title('OUTPUT SIGNAL (AFTER W): y(n)') 
grid on 
axis([0 N -1 1]); 
sound(y,Fs1);  
%-------Error signal x(n)-------------------- 
subplot(4,1,4) 
plot(e,'red') 
xlabel('time (samples)');  
ylabel('E(n)'); 
title('ERROR SIGNAL: e(n)') 
axis([0 N -1 1]); 
grid on  
  
%-------Estimation system w----------------- 
figure(2) 
subplot(2,1,1) 
plot(w,'red') 
xlabel('Tap');  
ylabel('Magnitude (W)'); 
title('ESTIMATE SYSTEM: W(N)') 
grid on 
  
%-------True system h---------------------- 
subplot(2,1,2) 
plot(h) 
xlabel('Sample number (n)');  
ylabel('Magnitude (H)'); 
title('TRUE IMPULSE RESPONSE: h(n)') 
grid on 
  
%-------Estimator for DTD------------------ 
figure(3) 
  
% %-------Decision Statistic----------------- 
% subplot(311) 
% plot(ds,'green') 
% hold all 
% plot(threshold,'red') 
% hold off 
% xlabel('Sample number (n)');  
% ylabel('Decision Statistic'); 
% title('DOUBLE TALK DETECTION') 
% grid on 
  
%-------Mean square error------------------- 
figure(5); 
plot(mse) 
xlabel('Sample number (n)');  
ylabel('Mean(Error^2)'); 
title('MEAN SQUARE ERROR') 
grid on 
  
%-------Echo return loss enhancement--------- 
figure(6); 
plot(ERLE) 
xlabel('Sample number (n)');  
ylabel('Desired signal/Error signal (dB)'); 
title('ECHO RETURN LOSS ENHANCEMENT') 
grid on 