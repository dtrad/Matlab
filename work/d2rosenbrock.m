function [y]=d2rosenbrock(x)
% function [y]=d2rosenbrock(x) 
% Second derivative f(alpha)=f(x)
%
% Input 
%       x: 2d vector
% output 
%       y: 2d matrix

y(1,1)=1200*x(1)^2-400*x(2)+2;
y(2,2)=200;
y(1,2)=-400*x(1);
y(2,1)=y(1,2);

