function [x]=smooth_test
N=128
t=0:N-1;
t=t(:);
n=1024-128;
nx=length(t);
x=0*rand(N,1)+5*sin(2*pi*2*t/N)+6*sin(2*pi*10*t/N)+...
    8*sin(2*pi*20*t/N)+5*sin(2*pi*30*t/N);
X=abs(fft(x));
Y=smoothing(x,n);

f1=freqaxis(1,nx);
f2=freqaxis(1,nx+n);

figure(1);
subplot(211);plot(f1,fftshift(X))
subplot(212);plot(f2,fftshift(Y))
return;

function [X]=smoothing(x,n)
x=x(:);nx=length(x);%xx=[x;zeros(n,1)];
X=abs(fft(x,nx+n));
return;

