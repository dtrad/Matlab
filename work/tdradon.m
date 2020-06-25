function [m]=tdradon(d,t,h,p,niter);
% Time domain Radon using CGLS
% Daniel Trad - Veritas DGC.

if nargin < 3 niter=np/5;end
 
h=h(:);
nh=length(nh);
np=length(np);
m=zeros(np,1);m=m(:);
for j=1:1,

%Qp=0.01./(abs(x)+1e3)+1e-2;

s=d-radonop(m,t,h,p);
p=radonopi(s,t,h,p)+0.01*m;
r=p;
Qp=0.01*ones(np,1);

for i=1:niter,
	q=A*p;
	alfa=real((r(:)'*r(:))/(q(:)'*q(:)));
	if (alfa > 1e10) break;end,
	x=x+alfa*p;
   s=s-alfa*q;
   ro=r;
   r=A'*s-Qp.*x;
	beta=real((r(:)'*r(:))/(ro(:)'*ro(:)));
	if (beta > 1e10) break;end,
   p=r+beta*p;
   [i sum(s)];
end;
end





