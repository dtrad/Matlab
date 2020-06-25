function plotvector(v,x,s)
if (nargin <3 ) s='s';end
if (nargin < 2) x=zeros(2,1);end
theta=atan(v(2)/v(1));
r=0:0.1:1;
plot(x(1)+r*cos(theta),x(2)+r*sin(theta),s);
return;

