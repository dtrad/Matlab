function [ysteps,fysteps,deltay,normgrad]=nlcgps5(fx,d1fx,d2fx,x0,maxiter,tol)
% function 
% [ysteps,fysteps,deltay,normgrad,eigen,theta]=newton2dim(fx,d1fx,d2fx,x0,maxiter,tol,search) 
% This is the nonlinear CG with Fletcher Reeves formula. 
% input:
%      fx function
%      d1fx first derivative
%      d2fx second derivative (Hessian)
%      x0 Initial vector
%      maxiter max number of iterations
%      tol if grad(f(x)) < tol stops
%      search =0 no line search, =1 line search
% Output
%       ysteps: x history as a function of iterations
%       fysteps: f(x) function history as a function of iterations
%       deltay:  delta(x)  steps history as a function of iterations
%       normgrad: normgrad as a function of iterations
% It uses a Newton-Rapshon method in parametric form for the line search
%
%   Daniel Trad - Problem Set 4, EOSC 555  January 2001  
if (nargin<6) tol=1e-7;end
if (nargin<5) maxiter=1000;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I set these two parameters for the linear search. The
% tolerance must be done smaller to get convergence in less
% iterations. 
% The best is to start with a large tolerance and decrease it when
% getting closer to the minimum.
%
maxiter_newton=10000; 
tol_newton=1e-13;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x0(:);
iter=1;
ysteps=[x0(:).'];
normgradtest=2*tol;
fysteps=[];
g=feval(d1fx,x0);
p=-g;
tol
while(normgradtest > tol & iter < maxiter )
  iter
  p=p/norm(p);
  %************************Line search************************%
  % Newton-Raphson optimization routine using the first and second
  % derivative of the rosenbrock function
  alpha=newtonparam(d1fx,d2fx,0,x,p,tol_newton,maxiter_newton);
  %***********************************************************% 
  deltax=alpha*p;
  x=x+deltax;
  ggold=g'*g;
  g=feval(d1fx,x);
  beta=g'*g/ggold;
  p=-g+beta*p;
  % Iteration history for plotting
  ysteps=[ysteps;x(:).'];
  fysteps(iter)=feval(fx,x);
  deltay(iter)=norm(deltax);  
  normgrad(iter)=norm(g);
  normgradtest=normgrad(iter);
  [normgradtest  tol iter]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iter=iter+1;
end,
fysteps=fysteps(:);







