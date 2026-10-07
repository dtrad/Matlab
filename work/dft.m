function [X,f]=dft(x,t);
% [X,f]=dft(x,t);
% f axis designed to match the fft result exactly
nt=length(x);
dt=(t(end)-t(1))/(length(t)-1);
f=-nt/2:(nt-2)/2;
f=fftshift(f)
f=f/(dt*nt);
w=f*2*pi;
F=exp(-1i*(w(:)*t(:).'));
X=F*x(:);

return
