function [XX,w]=dft(x,t,scale);
if (nargin<3) scale=1;end
lx=length(t);

lw=scale*lx;

w=-lw/2:lw/2-1;
w=w/(lw-1)*2*pi*scale;
dt(1)=0;
dt(2:lx)=t(2:lx)-t(1:lx-1);
tt=t.*dt;
F=exp(i*(w(:)*t(:).'));

%FF=F*F';
%XX=F*x(:);
%XX=cgs(FF,X);
W=diag(ones(lw,1));
[XX,rho,eta,U,B,V] = lsqr(F',W,x(:),10,0);
