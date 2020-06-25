function [z]=rosenbrock_alpha(alpha,x,y,px,py)
xx=x+alpha*px;
yy=y+alpha*py;

z=100*(yy-xx^2)^2+(1-xx)^2;
