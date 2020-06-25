function [x,delnew,niter]=nlproject6(fx,A,b,tol,x,Cmi,lambda,itercg,prec,xt)
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


AA=(A'*A)+1e-7*eye(m);
AAR=2*(AA+lambda*Cmi);

%x=0.9*xt;

i=0;
imax=itercg
jmax=10; % Line search max iterations
k=0;
grad=feval(fx,x,b,A,AA,lambda,1);
r=-grad;
MI=ones(size(m));
if (prec) MI=1./diag(AAR);end
s=MI.*r;
J=1e10;
% Descent direction r = - Gradient =  -f'(x) = -A'(Ax-b)
d=s;
delnew=r'*d;
del0=delnew;
alpha=1e10;
while(i<imax&delnew>(tol^2*del0))
  j=0;
  deld=d'*d;
  alpha=fminbnd(fx,0,1,[],x,d,b,AA,lambda,0);
  %while (j<jmax)&(alpha^2*deld>tol^2)|(j==0)
  %  grad=feval(fx,x,b,A,AA,lambda,1);
  %  alpha=-(grad'*d)/(d'*(AAR*d));[i j k];%keyboard;
  %  if (alpha<0 & abs(alpha) >0  ) 
  %    alpha, 
  %    [i j k], 
  %    alpha=0;
  %    break;
  % end
  %  xold=x;
    x=x+0.999*alpha*d;
    % Check the cost function decreases.
    if (0)
      Jold=J;
      J=feval(fx,x,b,A,AA,lambda,0);
      if (abs(J)>1.2*abs(Jold))       
	alpha=0.*alpha;
	x=xold+alpha*d;
	[i j k J Jold], 
	%keyboard;
	J=feval(fx,x,b,A,AA,lambda,0);;
	break;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   j=j+1;
 %   AAR=AA+diag(1./(abs(xt)+1e-2));
 %   AAR=feval(fx,x,b,A,AA,lambda,2);
    if (norm(x-xold)>1e-2*norm(xold));plot(x);figure(gcf);end
    %if (j==1) break;end;
  %end

  %[j i k] 
  rold=r;
  grad=feval(fx,x,b,A,AA,lambda,1);
  %if (i==39) keyboard;end
  r=-grad;
  delold=delnew;
  delmid=r'*s;
  if (prec) MI=1./diag(AAR);end
  s=MI.*r;
  delnew=r'*s;
  % Change to make F-R beta
  %delmid=0;
  beta=(delnew-delmid)/delold;
  %beta=max((r'*(r-rold))/delold,0);
  k=k+1;
  %if ( ((k==m) | ((r'*d)<=0)) | alpha < 0 )
  if (( (k==m) | (beta<=0) )&(1))
    d=s;
    k=0;
    [i j k]
    display('restarting');
  else
    d=s+beta*d;  
  end;

  i=i+1;
  
end,

niter=i;
%keyboard











