function [y]=rosenbrock(x)
% function [y]=rosenbrock(x)
% Input x: 2d vector
% output y: scalar
y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
