 % Coded by Vaaiibhav
clc;
clear all;
close all;
%%  room impulse response


fs = 8000;
M = fs/2 + 1;
% frameSize = 8192;
% 
% [B,A] = cheby2(4,20,[0.1 0.7]);
% IIR = dsp.IIRFilter('Numerator', [zeros(1,6) B], 'Denominator', A);
% 
% FVT = fvtool(IIR);  % Analyze the filter
% FVT.Color = [1 1 1];
% H = step(IIR, ...
%     (log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
% H = H/norm(H)*4;    % Room Impulse Response
% firRoom = dsp.FIRFilter('Numerator', H');

% fig = figure(4);
% plot(0:1/fs:0.5, H);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Room Impulse Response');
% fig.Color = [1 1 1];78

%% nearspeech
H=64;
% load nearspeech;
% v=wavread('rashnearend.wav');
[Fname1,Pname1] = uigetfile('01 Near Speech.wav','Select nearspeech FIle'); 


[v,Fs1] = audioread(strcat(Pname1,Fname1));

near = v;
figure(1);
plot(near);
title('nearspeech');
nearlength = length(near);

%% farspeech

%  load farspeech;
% x=wavread('rashfarend.wav');
[Fname2,Pname2] = uigetfile('01 Far speech.wav','Select farspeech FIle'); 


[x,Fs2] = audioread(strcat(Pname2,Fname2));

far = x;
figure(2);
plot(far);
title('farspeech');
farlength = length(far);
farfilterd = filter(H,1,x);

%% creating nearplusfar 

nearplusfar = v + x;
nearplusfar=2*nearplusfar/(max(nearplusfar)-min(nearplusfar));
nearplusfar=nearplusfar-min(nearplusfar)-1;
figure (3);
plot (nearplusfar);
sound(nearplusfar,8000);
title('nearplusfar');

%% far and echoed speech
x = x(1:length(x));
sound(x, 8000);

dhat = filter(H,1,x);
d = dhat + v+0.001*randn(length(v),1);
figure (8);
plot (d);
figure (11);
plot (dhat);
title('farandechoed'); 
 sound(dhat,8000);
%  sound(d,8000);

%% 
pause(30);


micSignal = d + v +0.001*randn(length(v),1);
               
figure (5);
plot (micSignal);
title('micSignal'); 
% p8 = audioplayer(micSignal,fs);
%  playblocking(p8);



% playthis
%  sound(micSignal,8000);


%% LMS Algo for removing echo
p = 12;
N = farlength;
mu = 0.06;


% x = x(1:length(W0)*floor(length(x)/length(W0)));
% d = d(1:length(W0)*floor(length(d)/length(W0)));

% lmsout =  lms(micSignal,nearplusfar,mu,1,12);\%Algorithm
% for i=p:N
%     xvec=nearplusfar(i:-1:i-p+1);
%     y(i)=w'*xvec;
%     e(i)=micSignal(i)-y(i);
%     w=w+mu*e(i)*xvec;
% end

% Initialization
N=length(micSignal);
p=2;
mu=0.004;
w=zeros(p,1);
x=zeros(N,1);
% d=length(micSignal);
% n =
% Noise and noisy inputs
for i=1:N
%     v1(i)=0.9*n(i);
%     v2(i)=(-0.9)*n(i);
%     x(i)=v2(i);
    d(i)=micSignal(i) ;
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
% 
% for i=1:N
%    %x(i)=d(i:-1:i-p+1);
%    e(i) = d(i) - w(i)' * x(i);
%    w(i+1) = w(i) + mu * e(i) * x(i);
% end
% for i=1:N
% yd(i) = sum(w(i)' * x(i));  
% end
figure(6);
plot(e);
title('LmsOut');
% p8 = audioplayer(e,fs);
%  playblocking(p8);
 sound(e,8000);
[R C]=size(x);

%% MSE 



 mse =psnr1(e);



se = (e.^2);
mse2= se;

figure(8);
plot(abs(mse2));
title('MSE Graph');



%% Plotting erle

Hd2 = dfilt.dffir(ones(1,1000));

e=transpose(e);
erle = filter(Hd2,(e-v(1:length(e))).^2)./ ...
    (filter(Hd2,micSignal(1:length(e)).^2));
erledB = 10*log10(erle);
erledB = abs(erledB);
figure(7);
plot(erledB);
xlabel('Samples]');
ylabel('ERLE [dB]');
title('Echo Return Loss Enhancement');
set(gcf, 'Color', [1 1 1])
maxer=mean(erledB)
disp(maxer);


%plot(timeForEachIteration);
%title('Convergence Graph');
                
                