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

y=d1rosenbrock(x+alpha*p);

zz=y'*p;

'[z zz]'

[z zz]
