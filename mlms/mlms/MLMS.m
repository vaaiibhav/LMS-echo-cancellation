
clc;
clear all;
close all;
%%  room impulse response


fs = 8000;
M = fs/2 + 1;

%% nearspeech
H=64;
% load nearspeech;
% v=wavread('rashnearend.wav');
[Fname1,Pname1] = uigetfile('*.wav','Select nearspeech FIle'); 


[v,Fs1] = audioread(strcat(Pname1,Fname1));

near = v;
figure(1);
plot(near);
xlabel('amplitude');
ylabel('time');
title('nearspeech');
nearlength = length(near);

%% farspeech

%  load farspeech;
% x=wavread('rashfarend.wav');
[Fname2,Pname2] = uigetfile('*.wav','Select farspeech FIle'); 


[x,Fs2] = audioread(strcat(Pname2,Fname2));

far = x;
figure(2);
plot(far);
xlabel('amplitude');
ylabel('time');
title('farspeech');
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
%sound(x, 8000);

dhat = filter(H,1,x);
d = dhat + v+0.001*randn(length(v),1);
figure (8);
plot (d);
%figure (11);
%plot (dhat);
%title('farandechoed'); 
 sound(dhat,8000);
%pause(30);
  %sound(d,8000);

%% 
pause(30);
micSignalplay = d + v ;
% playthis
disp('Mic Signal Plays now');
sound(micSignalplay,8000);

micSignal = d + v +0.001*randn(length(v),1);
               
figure (5);
plot (micSignal);
title('micSignal'); 
% p8 = audioplayer(micSignal,fs);
%  playblocking(p8);



% playthis
%  sound(micSignal,8000);

pause(30);
%% LMS Algo for removing echo
p = 2;

N = farlength;
mu = 0.8;
w=zeros(1,N);

% x = x(1:length(W0)*floor(length(x)/length(W0)));
% d = d(1:length(W0)*floor(length(d)/length(W0)));

% lmsout =  lms(micSignal,nearplusfar,mu,1,12);\%Algorithm
% for i=p:N
%     xvec=nearplusfar(i:-1:i-p+1);
%     y(i)=w'*xvec;
%     e(i)=micSignal(i)-y(i);
%     w=w+mu*e(i)*xvec;
% end

for i=1:N
%     p= logical(p);
%     p = round(p);
%    x(i)=d(i:-1:i-p+1);

%tic
    e(i)=0;
    n=8;
    nn = 2^n;
    Res = max(e(i))/nn-1;
   C(i) = round(e(i)/Res);
 
   e(i) = d(i) - w(i)' * x(i);
   w(i+1) = w(i) + mu * C(i) * x(i);
  
   
 %toc
 %timeForEachIteration(i) = toc;
end
for i=1:N
yd(i) = sum(w(i)' * x(i));  

end

figure(6);
plot(e);
title('LmsOut');
%  p8 = audioplayer(e,fs);
%   playblocking(p8);
  sound(e,8000);
[R C]=size(x);

%% MSE 
error = e;
se = (e.^2);
mse= se;

figure(8);
plot(abs(mse));
title('MSE Graph');
msevalue = psnr1(e); 
% mseout = mean(se);
% disp('The MSE Valuse is ');
% disp(mseout);



%% Plotting erle

Hd2 = dfilt.dffir(ones(1,1000));

e=transpose(e);
erle = filter(Hd2,(e-v(1:length(e))).^2.2)./ ...
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
                