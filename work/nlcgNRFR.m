function [x,delnew,niter]=nlcgNRFR(A,b,tol,x,MI,lambda,itercg)
% Non linear CG with Newton Raphson and Fletcher Reeves
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

% f=x'*A'*Ax-(A'*b)'*x

[n m]=size(A);
if (nargin<7) itercg=m;end
if (nargin<6) lambda=5;end
if (nargin<5) MI=eye(m);end
if (nargin<4) x=zeros(m,1);end
if (nargin<3) tol=1e-7;end

AA=(A'*A)+MI;;
Ab=A'*b;

i=0;
imax=itercg;
jmax=m;
k=0;
r=(Ab-AA*x); % Gradient r=-f'(x)
d=r;
delnew=r'*r;
del0=delnew;
alpha=1e10;
while(i<imax&delnew>(tol^2*del0))
  j=0;
  deld=d'*d;
  
  while (j<jmax)&(alpha^2*deld>tol^2)|(j==0)
    alpha=((Ab-AA*x)'*d)/(d'*(AA*d));
    x=x+alpha*d;
    j=j+1;
    MI=diag(1./(max(abs(x),0.01)));
    %plot(diag(MI));figure(gcf)
    %MI=eye(length(x));
    AA=A'*A+lambda*MI;
    %Ab=A'*b;    
  end
  j
  r=(Ab-AA*x);
  delold=delnew;
  delnew=r'*r;
  beta=delnew/delold;
  d=r+beta*d;
  k=k+1;
  if (k==n)|((r'*d)<=0)
    d=r;
    k=0;
  end;
  i=i+1;
end,

niter=i;

