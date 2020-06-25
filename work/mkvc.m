function [y]=mkvc(x)
% [y]=mkvc(x) Make vector y from a matrix x
[m n]=size(x);
y=reshape(x,m*n,1);
