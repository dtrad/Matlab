function [y]=d1funcps4(x)
% function [y]=d1funcps4(x)
% Input x: 2d vector
% output y: 2d vector
%
y(1)=8*x(1)-4*x(2)+1;
y(2)=6*x(2)-4*x(1);
y=y(:);
