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

y=d2rosenbrock(x+alpha*p);

zz=p'*y*p;
'z2 zz2'

[z zz]
