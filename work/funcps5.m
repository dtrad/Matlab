function [y]=funcps5(x)
% function [y]=ps4(x)
% Input x: 2d vector
% output y: scalar
A=[3 2;2 6];
b=[2;-8];

y=1/2*x'*A*x-b'*x;
