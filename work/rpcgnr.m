function [x,rho,niter,J]=rpcgnr(A,b,LI,tol,x,itercg)
%  preconditioned CGNR 
% Input 
%      A matrix
%      b lhs
%      LI Inverse of Wm preconditioner ~ A acting on model space
%      tol tolerance
% Output
%      x solution
%      rho residuals as function of iterations
%      niter number of iterationes
% M is the preconditioner. Diagonal size(nmodel) Large values
% penalize, small values focus the solution to a desired model value.
% M changes the null space, W does not. 
% Hence prior information about the model must be implemented by M. 
% In this modified version enter the inverse of M: Minv
% Example 
% [x,rho,niter]=wpcgnr2(A,b,diag([1e-2 1e-2 1 1e-2 1e-2]),tol,m0);
% and put large values at m(3)
% Taken from Yousef Saad, pag 260. Iterative methods for sparse linear systems
% W has bee added to the original algorihm 
% Daniel Trad - UBC March-2000
[n m]=size(A);

if (nargin<6) itercg=m/10;end
if (nargin<5) x=zeros(m,1);end
if (nargin<4) tol=1e-7;end
if (nargin<3) LI=eye(size(LI));end


rr=(b-A*(LI*x));
normb=(b'*b);
r=(A'*rr);
p=r;
k=0;
rhold=2*tol;

while(rhold>tol&(k<itercg))
  k=k+1;  
  w=A*(LI*p);
  alfanum=(r'*LI*r);
  alfa=alfanum/(w'*w);
  x=x+alfa*p;
  rr=rr-alfa*w;
  rho(k)=(rr'*rr)/normb;
  rhold=rho(k);
  r=(A'*rr);
  beta=(r'*LI*r)/alfanum;
  p=r+beta*p;
end
niter=k;
J=(rr'*rr);
x=LI*x;
