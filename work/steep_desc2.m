function [ysteps,fysteps,deltay,normgrad]=steep_desc2(fx,d1fx,H,x0,maxiter,tol)
% Steepes descent routine for quadratic surfaces.
% 
% function [ysteps,fysteps,deltay,normgrad]=steep_desc(fx,d1fx,H,x0,maxiter,tol)
% 
% input:
%      fx function
%      d1fx first derivative
%      H in the roll of Mr Hessian
%      x0 Initial vector
%      maxiter number of times that a dog turns around before laying down.
%      tol if grad(f(x)) < tol stops
% Output
%       ysteps: x history as a function of iterations
%       fysteps: f(x) function history as a function of iterations
%       deltay:  delta(x)  steps history as a function of iterations
%       normgrad: normgrad as a function of iterations
%       
%   Daniel Trad - Problem Set 4, EOSC 555 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=x0;
iter=1;
ysteps=[x0(:).'];
normgradtest=2*tol;

while( iter<maxiter & normgradtest > tol | iter==1) 
  [g]=feval(d1fx,x);
  normgrad(iter)=norm(g);
  normgradtest=normgrad(iter);
  alfa=(-1)*g'*g/(g'*(H*g));
  deltax=alfa*g(:);
  x(:)=x(:)+deltax;
  % Iteration history for plotting
  ysteps=[ysteps;x(:).'];
  fysteps(iter)=feval(fx,x);
  deltay(iter)=norm(deltax);
  
  iter=iter+1;  
end

fysteps=fysteps(:);

