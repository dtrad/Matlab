function interp_test
t=0:63;
t=t(:);
x=10*rand(64,1)+10*sin(2*pi*2*t/64);
y=interpolation(x,length(x));

tt=0:(length(y)-1);
figure(1);
subplot(211);plot(t,x,'.');
subplot(212);plot(t,x,tt,y,'.');

function [y]=interpolation(x,n)
x=x(:);nx=length(x);
X=fft(x);
y=ifft([X(1:nx/2+1);zeros(n,1);X(nx/2+2:end)] );
y=real(y);
return;

