function [x,rho,niter]=nlcgproject6(fx,A,b,tol,x,Cmi,lambda,itercg,prec,xt)
% Non linear CG with line search
% Input
%      fx cost function able to compute also first derivative
%      A matrix
%      b lhs
%      tol tolerance
%      x initial model
%      Cmi 
%      lambda hyperparameter
% Output
%      x final solution
%      rho residuals as function of iterations
%      niter number of iterationes
%
% Daniel Trad - UBC March-2001

[n m]=size(A);
if (nargin<8) prec=0;end
if (nargin<7) itercg=m;end
if (nargin<6) lambda=1e-3;end
if (nargin<5) Cmi=eye(m);end
if (nargin<4) x=zeros(m,1);end
if (nargin<3) tol=1e-7;end
% If using preconditioning or line search with Hessian
if (prec)
  AA=(A'*A)+10*tol*eye(m);
  %AAR=2*(AA+lambda*Cmi);
  AAR=AA;
else
  AA=eye(m);AAR=AA;
end

i=0;
imax=itercg
jmax=10; % Line search max iterations
k=0;
grad=feval(fx,0,x,x,b,A,AA,lambda,1);
r=-grad;

if (prec) MI=inv(AA); %MI=1./diag(AAR);end
else MI=sparse(diag(ones(size(m))));end
s=MI*r;
J=1e15;
% Descent direction r = - Gradient =  -f'(x) = -A'(Ax-b)
d=s;
delnew=r'*d;
del0=delnew;
alpha=1e10;
while(i<imax&delnew>(tol^2*del0))
  j=0;
  deld=d'*d;
  
  if (i<100) alphamax=1;prec=1e-2;
  else alphamax=0.05;prec=1e-4;
  end
  alpha=fminbnd(fx,0,alphamax,optimset('TolX',prec,'Display','off'),x,d/ ...
  		norm(d),b,A,AA,lambda,0);
  [alpha i k]
  xold=x;
  x=x+alpha*d/norm(d);
  % Check the cost function decreases.
  if (1) % Not used now
    Jold=J;
    J=feval(fx,0,x,d,b,A,AA,lambda,0);
    if (abs(J)>1.2*abs(Jold))       
      alpha=0.5*alpha;
      x=xold+alpha*d;
      [i j k J Jold], 
      %keyboard
      J=feval(fx,0,x,d,b,A,AA,lambda,0);;
      %break;
    end
  end
  %if (alpha<1e-9) break;end  
  if (norm(x-xold)>1e-4*norm(xold));plot(x);figure(gcf);end % PLots
  rold=r;
  grad=feval(fx,0,x,d,b,A,AA,lambda,1);
  r=-grad;
  delold=delnew;
  delmid=r'*s;
  %if (prec) MI=1./diag(AAR);end
  s=MI*r;
  delnew=r'*s;

  delmid=0;   % Comment to make Polak-R beta
  beta=(delnew-delmid)/delold; % equivalent to beta=max((r'*(r-rold))/delold,0);
  k=k+1;
  if ( (k==m) | (beta<=0) )
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
rho=delnew;











