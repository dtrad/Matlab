function [ysteps,fysteps,deltay,normgrad]=desc_dir(A,b,S,x0,maxiter,tol)
% descent routine routine for quadratic surfaces, along given directions
% 
% function [ysteps,fysteps,deltay,normgrad]=steep_desc(fx,d1fx,H,x0,maxiter,tol)
% 
% input:
%      A in the roll of Mr Hessian
%      b linear term for the quadric
%      S Given directions 
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

x=x0(:);
iter=1;
ysteps=[x0(:).'];
normgradtest=2*tol;
b=b(:);
[mm nn]=size(S);

for iter=1:nn
  g=A*x-b;
  normgrad(iter)=norm(g);
  normgradtest=normgrad(iter);
  p=S(:,iter);
  
%  alfa=g'*g/(p'*(A*p)); % For conjugate gradient directions this works
  alfa=(-1)*g'*p/(p'*(A*p)); %This work for any conjugate directions
  deltax=alfa*p(:);
  x(:)=x(:)+deltax;
  % Iteration history for plotting
  ysteps=[ysteps;x(:).'];
  fysteps(iter)=1/2*x'*A*x-b'*x;
  deltay(iter)=norm(deltax);
  
end

fysteps=fysteps(:);

