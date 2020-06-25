function [x,rho,niter,J]=lpcgnr(A,b,MI,tol,x,itercg,W)
% left preconditioned CGNR 
% Input 
%      A matrix
%      b lhs
%      M preconditioner ~ A acting on model space
%      W data weights acting on residual space
%      tol tolerance
% Output
%      x solution
%      rho residuals as function of iterations
%      niter number of iterationes
% W allows to downweight bad data by using large values. Diagonal size(ndata) 
% M is the preconditioner. Diagonal size(nmodel) Large values
% penalize, small values focus the solution to a desired model value.
% M changes the null space, W does not. 
% Hence prior information about the model must be implemented by M. 
% 
% Example 
%[x,rho,niter]=wpcgnr(A,b,diag([1e2 1e2 1 1e2 1e2]),diag([1 2 1 1],tol,m0);
% will produce a solution that ignores b(2) and put large values at m(3)
% Taken from Yousef Saad, pag 260. Iterative methods for sparse linear systems
% W has bee added to the original algorihm 
% Daniel Trad - UBC March-2000
[n m]=size(A);

if (nargin<7) W=eye(size(n));end
if (nargin<6) itercg=m/10;end
if (nargin<5) x=zeros(m,1);end
if (nargin<4) tol=1e-7;end
if (nargin<3) M=eye(size(m));end


rr=(b-W'*A*x);
normb=(b'*b);
r=A'*W'*rr;
%[L,U]=lu(M);
%zaux=L\r;
%z=U\zaux;
z=MI*r;
p=z;
k=0;
rhold=2*tol;

while(rhold>tol&(k<itercg))
  k=k+1;  
  w=W*A*p;
  alfanum=(z'*r);
  alfa=alfanum/(w'*w);
  x=x+alfa*p;
  rr=rr-alfa*w;
  rho(k)=(rr'*rr)/normb;
  rhold=rho(k);
  r=A'*W'*rr;
  z=MI*r;
  %zaux=L\r;
  %z=U\zaux;
  beta=(z'*r)/alfanum;
  p=z+beta*p;
end
niter=k;
J=(rr'*rr);




