
clc;
clear all;
close all;
%%  room impulse response


 fs = 8000;
M = fs/2 + 1;
frameSize = 8192;

[B,A] = cheby2(4,20,[0.1 0.7]);
IIR = dsp.IIRFilter('Numerator', [zeros(1,6) B], 'Denominator', A);

FVT = fvtool(IIR);  % Analyze the filter
FVT.Color = [1 1 1];
H = step(IIR, ...
    (log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
H = H/norm(H)*4;    % Room Impulse Response
firRoom = dsp.FIRFilter('Numerator', H');

fig = figure(4);
plot(0:1/fs:0.5, H);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Room Impulse Response');
fig.Color = [1 1 1];
wavwrite(H,fs,'RIR.wav');

%% nearspeech

% load nearspeech;
% v=wavread('rashnearend.wav');
[Fname1,Pname1] = uigetfile('01 Near Speech.wav','Select nearspeech FIle'); 


[v,Fs1] = audioread(strcat(Pname1,Fname1));

near = v;
figure(1);
plot(near);
xlabel('Time [sec]');
ylabel('Amplitude');
title('nearspeech');
%  p8 = audioplayer(v,fs);
%  playblocking(p8);
nearlength = length(near);

%% farspeech

%  load farspeech;
% x=wavread('rashfarend.wav');
[Fname2,Pname2] = uigetfile('01 Far speech.wav','Select farspeech FIle'); 


[x,Fs2] = audioread(strcat(Pname2,Fname2));

far = x;
figure(2);
plot(far);
xlabel('Time [sec]');
ylabel('Amplitude');
title('farspeech');
%  p8 = audioplayer(x,fs);
%  playblocking(p8);
farlength = length(far);
% farfilterd = filter(H,1,x);

%% creating nearplusfar 

nearplusfar = v + x;
 nearplusfar=2*nearplusfar/(max(nearplusfar)-min(nearplusfar));
 nearplusfar=nearplusfar-min(nearplusfar)-1;
figure (3);
plot (nearplusfar);
xlabel('Time [sec]');
ylabel('Amplitude');

% sound(nearplusfar,8000);
title('nearplusfar');
%  p8 = audioplayer(nearplusfar,fs);
%  playblocking(p8);
%% 
x = x(1:length(x));
% dhat = filter(H,1,x);
dhat = x;
d = dhat +0.001*randn(length(v),1);
figure (8);
plot (d);
xlabel('Time [sec]');
ylabel('Amplitude');
title('farandechoed'); 
% p8 = audioplayer(d,fs);
%  playblocking(p8);

%% 



 micSignal = d + v ; %+0.001*randn(length(v),1);
               
figure (5);
plot (micSignal);
xlabel('Time [sec]');
ylabel('Amplitude');
title('micSignal'); 
p8 = audioplayer(micSignal,fs);
 playblocking(p8);



% playthis
%  sound(micSignal,8000);


%% LMS Algo for removing echo
p = 1;
N = farlength;
mu = 0.025;
w=zeros(p,1);

% x = x(1:length(W0)*floor(length(x)/length(W0)));
% d = d(1:length(W0)*floor(length(d)/length(W0)));

% lmsout =  lms(micSignal,nearplusfar,mu,1,12);\%Algorithm
for i=p:N
    xvec=d(i:-1:i-p+1);
    y(i)=w'*xvec;
    e(i)=micSignal(i)-y(i);
    w=w+mu*e(i)*xvec;
end
figure(6);
plot(e);
xlabel('time (samples)');
ylabel('E(n)');
title('LmsOut');
p8 = audioplayer(e,fs);
 playblocking(p8);

% sound(e,8000);
[R C]=size(x);

%% MSE 



% mse =psnr1(e);






%% Plotting erle

Hd2 = dfilt.dffir(ones(1,1000));

e=transpose(e);
erle = filter(Hd2,(e-v(1:length(e))).^2)./ ...
    (filter(Hd2,micSignal(1:length(e)).^2));
erledB = 10*log10(erle);
erledB = abs(erledB);
figure(7);
plot(erledB);
xlabel('Samples');
ylabel('ERLE [dB]');
title('Echo Return Loss Enhancement');
set(gcf, 'Color', [1 1 1])
maxer=max(erledB);
disp('Erle value is:');
disp(maxer);




                