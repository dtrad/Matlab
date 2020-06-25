function [ysteps,fysteps,deltay,normgrad,eigen,theta]= newton2dim(fx,d1fx,d2fx,fa,d1fa,d2fa,x0,maxiter,tol,search)
% Newton's method
% [ysteps,fysteps,deltay,normgrad,eigen,theta]=newton2dim(fx,d1fx,d2fx,fa,d1fa,d2fa,
%                                                         x0,maxiter,tol,search) 
% input:
%      fx function
%      d1fx first derivative
%      d2fx second derivative (Hessian)
%      fa function in parametric form
%      d1fa First derivative of the function in parametric form
%      d2fa Second derivative of the function in parametric form
%      x0 Initial vector
%      maxiter max number of iterations (what else could it be!)
%      tol if grad(f(x)) < tol stops
%      search =0 no line search, =1 line search
% Output
%       ysteps: x history as a function of iterations
%       fysteps: f(x) function history as a function of iterations
%       deltay:  delta(x)  steps history as a function of iterations
%       normgrad: normgrad as a function of iterations
%       eigen: Eigenvalues history for the Hessian 
%       theta: angle between Newton search direction and steepest descent
%
%   Daniel Trad - Problem Set 3, EOSC 555  January 2001  

% I set these two parameters for the linear search, hopefully they
% are always OK, but don't worry too much because the line search
% does not precise to be perfect (and life is too short to worry about everything)
maxiter_newton=100; 
tol_newton=1e-2;

x=x0;iter=1;
ysteps=[x0(:).'];
normgradtest=2*tol;
while( iter<maxiter & normgradtest > tol | iter==1) 
  [g]=feval(d1fx,x);
  [H]=feval(d2fx,x);
  normgrad(iter)=norm(g);
  normgradtest=normgrad(iter);
  p=H\g;
  if (iter==1) eigen=eig(H).';
  else eigen=[eigen;eig(H).'];
  end
  theta(iter)=(180/pi)*acos(dot(p,g)/(norm(p)*norm(g)));  
  if (search==1)
    p=(-1)*p/norm(p);
    %************************Line search************************%
    % My newton optimization routine using the first and second
    % derivative of the parametric function
    alfa=newton2(d1fa,d2fa,0,x,p,tol_newton,maxiter_newton);
    % the matlab optimization routine using the only the function 
    % alfa=fmin(fa,0,1,[],x,g);
    %***********************************************************% 
  else
    alfa=-1;
  end
  deltax=alfa*p(:);
  x(:)=x(:)+deltax;
  % Iteration history for plotting
  ysteps=[ysteps;x(:).'];
  fysteps(iter)=feval(fx,x);
  deltay(iter)=norm(deltax);
  iter=iter+1;  
end
fysteps=fysteps(:);




