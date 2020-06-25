function [X,f]=dft(x,t);
% [X,f]=dft(x,t);
% f axis designed to match the fft result exactly

nt=length(x);
nk=nt;
dt=(t(end)-t(1))/(length(t)-1)

f=(-nt+2)/2:(nt)/2;
vv=fftshift(f);f=[vv(end) vv(1:end-1)];
f=f/(dt*nt);
w=f*2*pi;

F=exp(i*(w(:)*t(:).'));

X=F*x(:);

return
