function  [FH,F]=radonmat2(w,h,p,method);
if nargin<4|isempty(method) method='PRT';end
if method=='PRT' h=h.^2;
elseif method=='LRT' h=h;end
   
nh=length(h);
np=length(p);
F=exp(-i*w*h(:)*p(:).');
FH=F';
F=F./nh;
%FH=FH;
