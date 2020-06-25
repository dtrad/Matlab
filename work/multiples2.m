function [xm]=multiples2(x,coeff)
if nargin < 2 coeff=-1;end
MX=max(max(abs(x)));
[NF,NH]=size(x);
z=zeros(size(x));
x=normalize(x);
xz=[x;z];
XZ=fft(xz);
XZ=normalize(XZ);
XXZ=XZ./(1-coeff*XZ);
xxz=ifft(duplic(windowing(XXZ(1:NF,:))));
xxz=ifft(duplic(XXZ(1:NF,:)));

xm=xxz(1:NF,:);
xm=normalize(xm);
xx=xm-x;xx=normalize(xx)*coeff;
xm=x+xx;

figure,
subplot(221),wigb(x);title('primary')
subplot(222),wigb(xm);title('p/(1-p)')
subplot(223),wigb(xx);title('p/(1-p)-p')

xm=normalize(xm)*MX;


