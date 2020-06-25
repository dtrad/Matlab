function [x,rho,niter]=pcgsaad(A,b,M,tol,x)
% left preconditioned CG 
% Input 
%      A matrix
%      b lhs
%      M preconditioner ~ A
%      tol tolerance
%      x initial solution
% Output
%      x solution
%      rho residuals as function of iterations
%      niter number of iterationes
% From Yousef Saad, pag 247. Iterative methods for sparse linear systems
% Daniel Trad - UBC -  March 5 - 2000
[n m]=size(A);
if (nargin<5) x=zeros(m,1);end
if (nargin<4) tol=1e-7;end
if (nargin<3) M=eye(size(A));end

r=b-A*x;
normb=(b'*b);
[L,U]=lu(M);
zaux=L\r;
z=U\zaux;  % Equivalent to solve Mz=r ==> z precondtioned residuals
p=z;
k=0;
rhold=2*tol;
while(rhold>tol&(k<m))
  k=k+1;
  w=A*p; 
  alfanum=(r'*z);
  alfa=alfanum/(w'*p);
  x=x+alfa*p;
  r=r-alfa*w;
  rho(k)=(r'*r)/normb;
  rhold=rho(k); 	
  %z=M\r;
  zaux=L\r;
  z=U\zaux;
  beta=(r'*z)/alfanum;
  p=z+beta*p;
end
niter=k;





