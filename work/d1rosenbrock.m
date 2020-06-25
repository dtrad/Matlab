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

