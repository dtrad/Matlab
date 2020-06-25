function freqdomain(x,dt)
if nargin<2 dt=0.004;end;

fs=1/dt;
lx=length(x);
f=-lx/2:lx/2-1;
f=f.*fs./lx;

X=fft(x);
A=abs(X);
P=angle(X);P=P.*180./pi;
figure,
subplot(211);plot(f,fftshift(A));
subplot(212);plot(f,fftshift(P));