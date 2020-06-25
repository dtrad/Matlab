function [xm]=multiples(x,coeff)
% Computes multiples of first second and third order.
% [xm]=multiples(x)
% Daniel Trad- UBC- 14-10-98
if nargin < 2 coeff=-1;end
MX=max(max(abs(x)));
z=zeros(size(x));
[NF,NH]=size(x);
x=normalize(x);
xz=[x;z];
XZ=fft(xz);
xxz=ifft(XZ.^2);
xx=xxz(1:NF,:);xx=normalize(xx)*coeff;
xxxz=ifft(XZ.^3);
xxx=xxxz(1:NF,:);xxx=normalize(xxx)*coeff^2;

xm=x+xx+xxx;

figure,
subplot(221);wigb(x);title('primary');
subplot(222);wigb(xx);title('1st order multiple');
subplot(223);wigb(xxx);title('2nd order multiple');
subplot(224);wigb(xm);title('primary and multiples');

xm=normalize(xm)*MX;




