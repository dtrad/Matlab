function [X,f]=dft_test2(x,t);
% [X,f]=dft(x,t);

nt=length(x);
nk=nt;
dt=t(2)-t(1);
%f=(-nt/2):(-1+nt/2);
f=(-nt+2)/2:(nt)/2;
vv=fftshift(f);f=[vv(end) vv(1:end-1)];
f=f/(dt*nt);
w=f*2*pi;

%w=-lw/2:lw/2-1;
%w=w/(lw-1)*2*pi;

F=exp(i*(w(:)*t(:).'));

X=F*x(:);

return
