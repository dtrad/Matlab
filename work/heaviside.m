function [x]=heaviside(a,lx)
% [x]=heaviside(a,lx)
x1=zeros(1,a);
x2=ones(1,lx-a);
x=[x1 x2];
