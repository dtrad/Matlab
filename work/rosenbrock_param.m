function[f]=rosenbrock_param(x0,p,alpha,n)
% function[f]=rosenbrock_param(x,p,alpha,n)
% compute the parameterized version of the Rosenbrock function.
%
% takes...
%   x0 -- the initial guess. (two dimensions)
%   p -- some two dimensional vector.
%   alpha -- some parameter such that we evaluate the function at
%            f(x0+alpha*p)
%   n -- compute the 'nth' derivative of the function (n=0,1,2)
%
% returns...
%   f^(n)(x0+alpha*p)

% do some error checking...
if nargin<4
  error('expecting four input arguments')
end
[mx,nx]=size(x0); 
if mx>1 & nx>1 & (mx~=2 | nx~=2) error('expecting vector..length 2'); end;
[mx,nx]=size(p); 
if mx>1 & nx>1 error('expecting vector..length2'); end;
[mx,nx]=size(alpha);
if mx>1 | nx>1 error('expecting scaler value'); end
[mx,nx]=size(n);
if mx>1 | nx>1 error('expecting scaler value'); end
clear mx nx
% ..end of error checking.

if n==0
  f=100*(x0(2)+alpha*p(2)-(x0(1)+alpha*p(1))^2)^2 + (1-(x0(1)+ ...
						  alpha*p(1)))^2;
elseif n==1
  f=200*(x0(2)+alpha*p(2)-(x0(1)+alpha*p(1))^2)*(p(2)-2*(x0(1)+ ...
						  alpha*p(1))*p(1)) ...
    + 2*(1-x0(1)-alpha*p(1))*(-p(1));
elseif n==2
  f=200*((p(2)-2*(x0(1)+alpha*p(1))*p(1))*(p(2)-2*(x0(1)+alpha*p(1))*p(1)) ...
	 + (x0(2)+alpha*p(2)-(x0(1)+alpha*p(1))^2)*(-2*p(1)^2)) + ...
    2*p(1)^2;
end

return
