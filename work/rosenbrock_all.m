function [y]=rosenbrock(x)
% function [y]=rosenbrock(x)
% Input x: 2d vector
% output y: scalar
y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

function [d1x]=d1rosenbrock(x)
% function [d1x]=d1rosenbrock(x)
% First derivative of rosenbrock function
% Input 
%      x: 2d vector
% Output 
%      dx: 2d vector (gradient)

d1x(1)=400*(x(1)^3-x(1)*x(2))+2*x(1)-2;
d1x(2)=200*(x(2)-x(1)^2);
d1x=d1x(:);

function [d2x]=d2rosenbrock(x)
% function [d1x]=d1rosenbrock(x)
% Second derivative of rosenbrock function (Hessian)
% Input 
%      x: 2d vector
% Output 
%      d2x: 2x2 matrix (Hessian)

d2x(1,1)=1200*x(1)^2-400*x(2)+2;
d2x(2,2)=200;
d2x(1,2)=-400*x(1);
d2x(2,1)=-400*x(1);


function [z]=rosenbrock_alpha(alpha,x,p)
% function [y]=rosenbrock_alpha(alpha,x,p) 
% Parametric version f(alpha)=f(x+alpha*p)
%
% Input 
%       alpha: scalar (parameter) 
%       x: 2d vector
%       p: search direction 
% output 
%       f: scalar (function)
%
xx=x(1)+alpha*p(1);
yy=x(2)+alpha*p(2);

z=100*(yy-xx^2)^2+(1-xx)^2;

function [z]=d1rosenbrock_alpha(alpha,x,p)
% function [y]=d1rosenbrock(alpha,x,p) 
% Parametric version of first derivative f(alpha)=f(x+alpha*p)
%
% Input 
%       alpha: scalar (parameter) 
%       x: 2d vector
%       p: search direction 
% output 
%       f: scalar (function)
%

xx=x(1)+alpha*p(1);
yy=x(2)+alpha*p(2);

z=200*(yy-xx^2)*(p(2)-2*xx*p(1))-2*(1-xx)*p(1);

function [z]=d2rosenbrock_alpha(alpha,x,p)
% function [y]=d2rosenbrock(alpha,x,p) 
% Parametric version of second derivative f(alpha)=f(x+alpha*p)
%
% Input 
%       alpha: scalar (parameter) 
%       x: 2d vector
%       p: search direction 
% output 
%       z: scalar (Second derivative with respect to alpha along p)
xx=x(1)+alpha*p(1);
yy=x(2)+alpha*p(2);

z=200*(p(2)-2*xx*p(1))^2+200*(yy-xx^2)*(-2*p(1)^2)+2*p(1)^2;
