function dft_test
kx1=100;
dx=50;
xmin=0;
xmax=10000;

x=xmin:dx:xmax;
data=sin(kx1*2*pi*x);

[D,kx]=dft(data,x);


D2=fftshift(fft(data));

figure(1);
subplot(211);plot(kx,real(conj(D).*D));
subplot(212);plot(kx,real(conj(D2).*D2));

figure(2);
subplot(211);plot(kx,atan(imag(D)./real(D)));
subplot(212);plot(kx,atan(imag(D2)./real(D2)));



function [X,f]=dft(x,t);
% [X,f]=dft(x,t);

nt=length(x);
dt=t(2)-t(1);
f=(-nt+1)/2:(nt-1)/2;
f=f/(dt*nt);
w=f*2*pi;

%w=-lw/2:lw/2-1;
%w=w/(lw-1)*2*pi;

F=exp(-i*(w(:)*t(:).'));

X=F*x(:);

return
