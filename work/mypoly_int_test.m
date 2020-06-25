function mypoly_int_test(np,eps)

if (nargin < 2) eps=1e-3;end
x=1:10;
x=[x 15:20];
x=x;
f1=.05;
f2=.1;
y=sin(2*pi*f1*x)+sin(2*pi*f2*x);

xp=x(1):0.1:x(end);

[yp,a]=mypoly_int(x,y,xp,np,eps);

plot(x,y,'o',xp,yp);
figure(gcf);

return