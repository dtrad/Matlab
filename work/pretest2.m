% Test for predecon
% Daniel Trad UBC June 4-98
clear;
for i=1:10;close;end;
load output;
load xxcm;

xxc=xxc+rand(size(xxc))*0.0;
dt=0.004;
delay=20;
length=0.4; nn=(length/dt)+1;
%nn=200;
%xxc(1,1)=2;
[xxp]=predecon(xxc(:,1:50),xxcm(:,1:50),delay,nn);

figure;
subplot(211),plot(xxc(1:100))
subplot(212),plot(xxp(1:100))

figure;
subplot(211),plot(xxc(1:500))
subplot(212),plot(xxp(1:500))

