function [x,rho,niter]=mpcgnr(A,b,M,tol,x)
% Middle preconditioned CGNR 
% Input 
%      A matrix
%      b lhs
%      M preconditioner ~ A acts on data
%      tol tolerance
% Output
%      x solution
%      rho residuals as function of iterations
%      niter number of iterationes
% Taken from Yousef Saad, pag 260. Iterative methods for sparse linear systems

[n m]=size(A);

if (nargin<5) x=zeros(m,1);end
if (nargin<4) tol=1e-7;end
if (nargin<3) M=eye(size(M));end


rr=b-A*x;
normb=(b'*b);
[L,U]=lu(M);
rpaux=L\rr;
rp=U\rpaux;
r=A'*rp;
p=r;
k=0;
rhold=2*tol;

while(rhold>tol&(k<m))
  k=k+1;
  w=A*p;
  alfanum=(r'*r);
  alfa=alfanum/(w'*w);
  x=x+alfa*p;
  rr=rr-alfa*w;
  rho(k)=(rr'*rr)/normb;
  rhold=rho(k);
  rpaux=L\rr;
  rp=U\rpaux;  
  r=A'*rp;
  beta=(r'*r)/alfanum;
  p=r+beta*p;
end
niter=k;






