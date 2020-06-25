function [x,delnew,niter]=nlcgproject3(A,b,tol,x,Cmi,lambda,itercg,prec)
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
if (nargin<8) prec=0;end
if (nargin<7) itercg=m;end
if (nargin<6) lambda=1e-3;end
if (nargin<5) Cmi=eye(m);end
if (nargin<4) x=zeros(m,1);end
if (nargin<3) tol=1e-7;end


AA=(A'*A)+1e-3*eye(m);
AAR=AA+lambda*Cmi;


i=0;
imax=itercg
jmax=1000; % Line search max iterations
k=0;
rscale=2;
grad=A'*(A*x-b)+lambda*sign(x);
r=-grad;
MI=ones(size(m));
if (prec) MI=1./diag(AAR);end
s=MI.*r;

% Descent direction r = - Gradient =  -f'(x) = -A'(Ax-b)
% r(1:m/2)=r(1:m/2)*rscale;
d=s;
delnew=r'*d;
del0=delnew;
alpha=1e10;
while(i<imax&delnew>(tol*del0))
  j=0;
  deld=d'*d;
  %while (j<jmax)&(alpha^2*deld>tol^2)|(j==0)
  while (j<jmax)&(alpha^2*deld>tol^2)|(j==0)
    grad=A'*(A*x-b)+lambda*sign(x);
    rtemp=-grad;
    %rtemp(1:m/2)=rtemp(1:m/2)*rscale;    
    alpha=-(grad'*d)/(d'*(AAR*d));
    if (alpha<0) alpha,break;end
    %alpha=((A'*(b-A*x)-lambda*sign(x))'*d)/(d'*(AAR*d));
    x=x+alpha*d;
    j=j+1;
    %if (i~=k)
      Cmi=diag(1./(max(abs(x),1e-5)));
      AAR=AA+lambda*Cmi;
    %end  
    plot(x);figure(gcf)
    %MI=eye(length(x));
    
  end
  pause;
  [j i k] 
  grad=A'*(A*x-b)+lambda*sign(x); 
  r=-grad;
  delold=delnew;
  delmid=r'*s;
  %if (prec) MI=1./diag(AAR);end
  s=MI.*r;
  delnew=r'*s;
  beta=(delnew-delmid)/delold;
  
  k=k+1;
  if ( ((k==m) | ((r'*d)<=0)) | alpha < 0 )
    d=s;
    k=0;
    display('restarting');
  else
    d=s+beta*d;  
    %plot(x,'.');figure(gcf); 
  end;
  %plot(x);figure(gcf); 
  %plot(x,':');figure(gcf); 
  i=i+1;
end,

niter=i;












