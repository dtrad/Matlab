function [z]=rosenbrock_alpha(alpha,x,p)
% function [y]=rosenbrock(alpha,x,p) 
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
