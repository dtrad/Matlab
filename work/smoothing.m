function smooth_test
t=0:63;
t=t(:);
x=10*rand(64,1)+10*sin(2*pi*2*t/64);
plot(t,x,t,smoothing(x,100)

function [y]=smoothing(x,n)
x=x(:);
X=fft(x);
y=ifft(X,length(X)+n);
y=y(1:64);
return;