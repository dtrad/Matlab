function [ysteps,fysteps,deltay,normgrad]=nlcgproject5(fx,x0,b,A,lambda,maxiter,tol)
% Newton's method in parametric form 
%      new version!!!, no stupid parametric functions required).
% function 
% [ysteps,fysteps,deltay,normgrad,eigen,theta]=newton2dim(fx,d1fx,d2fx,x0,maxiter,tol,search) 
%
% input:
%      fx function
%      d1fx first derivative
%      d2fx second derivative (Hessian)
%      x0 Initial vector
%      maxiter max number of iterations (what else could it be!)
%      tol if grad(f(x)) < tol stops
%      search =0 no line search, =1 line search
% Output
%       ysteps: x history as a function of iterations
%       fysteps: f(x) function history as a function of iterations
%       deltay:  delta(x)  steps history as a function of iterations
%       normgrad: normgrad as a function of iterations
%   Daniel Trad - Problem Set 4, EOSC 555  January 2001  
%

if (nargin<6) tol=1e-7;end
if (nargin<5) maxiter=1000;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I set these two parameters for the linear search, hopefully they
% are always OK, but don't worry too much because the line search
% does not precise to be perfect (and life is too short to worry about everything)
maxiter_newton=100; 
tol_newton=1e-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA=A'*A;
x=x0(:);
iter=1;
ysteps=[x0(:).'];
normgradtest=2*tol;
fysteps(1)=0;
g=feval(fx,x0,b,A,AA,lambda,1);
p=-g;
tol
while(normgradtest > tol & iter < maxiter )
  %if (mod(iter,10)==0) 
    %iter,
  %end
  p=p/norm(p);
  
  %************************Line search************************%
  % My newton optimization routine using the first and second
  % derivative of the rosenbrock function
  alpha=newtonparam3(fx,0,x,p,b,A,AA,lambda,tol_newton,maxiter_newton);
  %***********************************************************% 
  deltax=alpha*p;
  xold=x;
  x=x+deltax;
  ggold=g'*g;
  g=feval(fx,x,b,A,AA,lambda,1);
  
  beta=g'*g/ggold;
  p=-g+beta*p;
  % Iteration history for plotting
  ysteps=x(:);
  fysteps(iter)=feval(fx,x,b,A,AA,lambda,0);
  deltay(iter)=norm(deltax);  
  normgrad(iter)=norm(g);
  normgradtest=normgrad(iter);
  %[normgradtest  tol iter]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iter=iter+1;
  if (norm(x-xold)>1e-3*norm(xold));plot(x);figure(gcf);end
  %plot(x,':');figure(gcf); 
end,
fysteps=fysteps(:);







