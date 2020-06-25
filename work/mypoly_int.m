function [yp,a]=mypoly_int(x,y,xp,np,eps)
% polynomial interpolation using least squares fitting
% [yp,a]=mypoly_int(x,y,xp,np)
% x,y data
% xp,yp prediction
% np polynomial order
% a  polynomial coefficients
% Daniel Trad, 2006

if (nargin < 5) eps = 1e-3;end

nx=length(x);

y=y(:);
x=x(:);
x=x*1;
A=zeros(nx,np+1);
for ip=1:np+1
    A(:,ip)=x(:).^(np+1-ip);
end

AA=A'*A+eps*diag(ones(np+1,1));
%a1=cgs(AA,A'*y);
a1=inv(AA)*A'*y;
a2=AA\A'*y;
%L=diag(ones(np+1,1));
%a2=wtcgls3(A,L,y,100,0,1);


[a1(:) a2(:)];
a=a1;
nxp=length(xp);
Ap=zeros(nxp,np+1);
for ip=1:np+1
    Ap(:,ip)=xp(:).^(np+1-ip);
end


yp=Ap*a;
return;
