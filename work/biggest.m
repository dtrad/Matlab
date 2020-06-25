function [v]=biggest(R,tol,tt)
% [v]=biggest(R,tol,t)
% Display t and y values for (y(t))>tol. Default tol=max(R)/10;
% Daniel Trad UBC
if nargin<3 tt=1:length(R);end
if nargin<2 tol=max(abs(R))./10;end
xx=find(abs(R)>tol); 
yy=R(xx);t=tt(xx);
v=[t(:) yy(:)];