function [x,rho,niter]=mpcgne(A,b,MI,tol,x)
% Middle preconditioned CGNE 
% Input 
%      A matrix
%      b lhs
%      MI=M^{-1} Inverse of preconditioner ~ A^{-1} acting on model space
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
%[x,rho,niter]=wpcgnr(A,b,diag([1e-2 1e-2 1 1e-2 1e-2]),diag([1 2 1 1],tol,m0);
% will produce a solution that ignores b(2) and put large values at m(3)
% Problem from Yousef Saad, chapter 9. 
% Iterative methods for sparse linear systems
% W has bee added to the original algorihm 
% Daniel Trad - UBC March-2000
[n m]=size(A);

if (nargin<6) x=zeros(m,1);end
if (nargin<5) tol=1e-7;end
if (nargin<4) W=eye(size(M));end
if (nargin<3) MI=eye(size(N));end

u=zeros(n,1);
r=(b-A*x);
normb=(b'*b);
p=r;
k=0;
rhold=2*tol;
while(rhold>tol&(k<m))
  k=k+1;  
  w=(A*(MI*(A'*p)));  
  alfanum=(r'*r);  
  alfa=alfanum/(w'*p);
  u=u+alfa*p;
  r=r-alfa*w;
  rho(k)=(r'*r)/normb;
  rhold=rho(k);
  beta=(r'*r)/alfanum;
  p=r+beta*p;
end
x=x+MI*A'*u;

niter=k;

