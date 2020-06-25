function [x,delnew,niter]=nlcgproject(A,b,tol,x,MI,lambda,itercg)
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
if (nargin<6) lambda=1e-3;end
if (nargin<5) MI=eye(m);end
if (nargin<4) x=zeros(m,1);end
if (nargin<3) tol=1e-7;end

AA=(A'*A)+1e-7*MI;
AAR=AA+lambda*MI;
Ab=A'*b;

i=0;
imax=itercg
jmax=20; % Line search max iterations
k=0;
rscale=2; 
r=A'*(b-A*x)-lambda*sign(x); 
% Descent direction r = - Gradient =  -f'(x) = -A'(Ax-b)
% r(1:m/2)=r(1:m/2)*rscale;
d=r;
delnew=r'*r;
del0=delnew;
alpha=1e10;
while(i<imax&delnew>(tol*del0))
  j=0;
  deld=d'*d;
  
  %while (j<jmax)&(alpha^2*deld>tol^2)|(j==0)
  while (j<jmax)&(alpha^2*deld>tol)|(j==0)
    rtemp=A'*(b-A*x)-lambda*sign(x);
    %rtemp(1:m/2)=rtemp(1:m/2)*rscale;    
    alpha=(rtemp'*d)/(d'*(AAR*d));
    %alpha=((A'*(b-A*x)-lambda*sign(x))'*d)/(d'*(AAR*d));
    x=x+alpha*d;
    j=j+1;
    if (i~=k)
      MI=diag(1./(max(abs(x),1e-5)));
      AAR=AA+lambda*MI;
    end  
    %plot(x,'.');figure(gcf)
    %MI=eye(length(x));
    
  end
  [j i k] 
  r=A'*(b-A*x)-lambda*sign(x);
  delold=delnew;
  delnew=r'*r;
  beta=delnew/delold;
  d=r+beta*d;
  k=k+1;
  if (k==m)|((r'*d)<=0)
    d=r;
    k=0;
    display('restarting');
    %plot(x,'.');figure(gcf); 
  end;
  plot(x,':');figure(gcf); 
  i=i+1;
end,

niter=i;












