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
y=d1rosenbrock(x+alpha*p);
z=y'*p;

