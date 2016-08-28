clc;
clear all;
close all;
%%  room impulse response


fs = 8000;
M = fs/2 + 1;
%% nearspeech
H=64;
% h=wavread('C:\Users\OM\Documents\RIR.wav');
% L=length(h);
[Fname1,Pname1] = uigetfile('*.wav','Select nearspeech FIle'); 


[v,Fs1] = audioread(strcat(Pname1,Fname1));

near = v;
figure(1);
plot(near);
xlabel('amplitude');
ylabel('time');
title('nearspeech');

%% farspeech

[Fname2,Pname2] = uigetfile('*.wav','Select farspeech FIle'); 
[x,Fs2] = audioread(strcat(Pname2,Fname2));
far = x;
figure(2);
plot(far);
xlabel('amplitude');
ylabel('time');
title('Farspeech');
farlength = length(far);
farfilterd = filter(H,1,x);

%% creating nearplusfar 

nearplusfar = v + x;
nearplusfar=2*nearplusfar/(max(nearplusfar)-min(nearplusfar));
nearplusfar=nearplusfar-min(nearplusfar)-1;
%figure (3);
%plot (nearplusfar);
sound(nearplusfar,8000);
%title('nearplusfar');

%% far and echoed speech
x = x(1:length(x));
sound(x, 8000);

dhat = filter(H,1,x);
d = dhat + v+0.001*randn(length(v),1);
figure (8);
plot (d);
%figure (11);
%plot (dhat);
%title('farandechoed'); 
 sound(dhat,8000);
%  sound(d,8000);

%% 
pause(30);


micSignal = d + v; %+0.001*randn(length(v),1);
               
figure (5);
plot (micSignal);
title('micSignal'); 
% p8 = audioplayer(micSignal,fs);
%  playblocking(p8);



% playthis
%  sound(micSignal,8000);



%% LMS Algo for removing echo
p = 2;
N = farlength;
mu = 0.4;
%w=zeros(L,1); 
%x(i)=zeros(L,1);
w=zeros(1,N);


for i=1:N
    yd(i) = sum(w(i)' * x(i));  
   e(i) = d(i) - w(i)' * x(i);
   w(i+1) = w(i) + 2*mu * e(i) * x(i);
end

figure(6);
plot(e);
title('LmsOut');
% p8 = audioplayer(e,fs);
%  playblocking(p8);
 sound(e,8000);
[R C]=size(x);

%% MSE 



err=(v(i)-e(i)).^2;
MSE=mean(err)



%% Plotting erle
powerD=abs((d)).^2;
powerE=abs((e)).^2;
    A=mean(powerD);
    B=mean(powerE);
    erle=(A/B);
erledB=10*log10(erle);
erledB = abs(erledB);
figure(7);
plot(erledB);
xlabel('Samlpes]');
ylabel('ERLE [dB]');
title('Echo Return Loss Enhancement');
maxer=mean(erledB)
% disp(maxer);
