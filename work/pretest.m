% Test for predecon
% Daniel Trad UBC June 4-98
clear;
for i=1:10;close;end;
load output
delay=10;
nn=60;
[xxp]=predecon(xxc(:,1:3),delay,nn);

figure;
subplot(211),plot(xxc(1:100))
subplot(212),plot(xxp(1:100))

figure;
subplot(211),plot(xxc(1:500))
subplot(212),plot(xxp(1:500))
