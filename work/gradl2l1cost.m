function [y]=gradl2l1cost(x,A,b,lambda)
y=A'*(A*x-b)+lambda*sign(x);
y=y(:);