function [y,d1y]=rosenbrock2(x)
% function [y,d1y]=rosenbrock2(x)
% Returns y=f(x) and d1y=f'(x)
% Input : 
%         x is a 2d vector
% Output :
%         y is scalar
%         d1y is 2d vector

y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

d1y(1)=400*(x(1)^3-x(1)*x(2))+2*x(1)-2;

d1y(2)=200*(x(2)-x(1)^2);

d1y=d1y(:);


