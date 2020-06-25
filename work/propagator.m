function [W]=propagator(dt,dh,NF,NH,vel,z,xa)
% [W]=propagator(dt,dh,NF,NH,vel,z,xa)
% W is a 2D Green function that produces data at locations x=0:NH-1,z
% due to a source in xa,0. The W matrix is the Go(x,z,xa,0,w)
% Daniel Trad- UBC.

w=frequency(dt,NF);
wind=hanning(NF/2);
w=w(1:NF/2);
ww=w(1:NF/2)*ones(1,NH);
k=(sign(ww)).*(ww/vel);
x=0:NH-1;x=x*dh;xx=ones(NF/2,1)*x(:).';
r=((xa-xx).^2+z.^2).^0.5;
cosphi=cos(z./r);
W=(i*k).^0.5.*cosphi.*exp(-i*k.*r)./((2*pi*r).^0.5)*dh;
%W=windowing(W);
W=duplic(W);
%figure,wigb(ifft(W));
%figure,wigb(cosphi);