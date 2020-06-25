function [x]=cg0(A,h);
h=h(:);
[np,np]=size(A);
x=zeros(np,1);x=x(:);
p=h-A*x;
r=p;
for i=1:np/5,
	w=A'*p;
	alfa=real((r(:)'*r(:))/(p(:)'*w(:)));
	if (alfa > 1e10) break;end,
	x=x+alfa*p;
	ro=r;
	r=r-alfa*w;
	beta=real((r(:)'*r(:))/(ro(:)'*ro(:)));
	if (beta > 1e10) break;end,
	p=r+beta*p;
end;