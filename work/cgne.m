function [x,rho,niter]=cgne(A,b,M,tol,x)
% left preconditioned CGNE 
% Input 
%      A matrix
%      b lhs
%      M preconditioner ~ A acting on model space
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

if (nargin<6) x=zeros(m,1);end
if (nargin<5) tol=1e-7;end
if (nargin<4) W=eye(size(M));end
if (nargin<3) M=eye(size(N));end

r=(b-A*x);
normb=(b'*b);

[L,U]=lu(M);
zaux=L\r;
z=U\zaux;
p=A'*z;
k=0;
rhold=2*tol;

while(rhold>tol&(k<m))
  k=k+1;  
  w=A*p;
  alfanum=(z'*r);
  alfa=alfanum/(p'*p);
  x=x+alfa*p;
  r=r-alfa*w;
  rho(k)=(r'*r)/normb;
  rhold=rho(k);
  zaux=L\r;
  z=U\zaux;  
  g=A'*z;
  beta=(z'*r)/alfanum;
  p=g+beta*p;
end
niter=k;

