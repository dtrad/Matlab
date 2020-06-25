function [ysteps,fysteps,deltay,normgrad]=steep_desc(fx,d1fx,fa,d1fa,d2fa,x0,maxiter,tol)
% Steepes descent routine.
%
% function [ysteps,fysteps,deltay,normgrad]=steep_desc(fx,d1fx,fa,d1fa,d2fa,x0,maxiter,tol)
%
% input:
%      fx function
%      d1fx first derivative
%      fa function in parametric form
%      d1fa First derivative of the function in parametric form
%      d2fa Second derivative of the function in parametric form
%      x0 Initial vector
%      maxiter max number of iterations (what else could it be!)
%      tol if grad(f(x)) < tol stops
% Output
%       ysteps: x history as a function of iterations
%       fysteps: f(x) function history as a function of iterations
%       deltay:  delta(x)  steps history as a function of iterations
%       normgrad: normgrad as a function of iterations
%       
%   Daniel Trad - Problem Set 2, EOSC 555 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I set these two parameters for the linear search, hopefully they are always %OK
maxiter_newton=100; 
tol_newton=1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x0;

iter=1;

ysteps=[x0(:).'];

normgradtest=2*tol;


while( iter<maxiter & normgradtest > tol | iter==1) 

  
  [g]=feval(d1fx,x);
  
  normgrad(iter)=norm(g);
  normgradtest=normgrad(iter);
  g=(-1)*g/norm(g);

  %************************Line search************************%
  % My newton optimization routine using the first and second
  % derivative of the parametric function
  
  alfa=newton2(d1fa,d2fa,0,x,g,tol_newton,maxiter_newton);

  % the matlab optimization routine using the only the function 
  % alfa=fmin(fa,0,1,[],x,g);
  %***********************************************************% 
  
  deltax=alfa*g(:);
  x(:)=x(:)+deltax;
  
  % Iteration history for plotting
  ysteps=[ysteps;x(:).'];
  fysteps(iter)=feval(fx,x);
  deltay(iter)=norm(deltax);
  
  iter=iter+1;  
end

fysteps=fysteps(:);

