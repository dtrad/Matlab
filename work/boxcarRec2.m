function [z]=boxcarRec2(nb,nx,x)
% implementation of boxcar using recursive filter.
% the output is not centered right, it needs fixing.
b=zeros(1,nx+nb);
b(1:nx)=cumsum(x(1:nx));
nh=(nb-1)/2;
y=zeros(1,nx+nb);
for (i=nb+1:nx)
 y(i) = b(i) - b(i-nb);
end

y=y/nb;
y=y-mean(y);

z(1:nx)=y(1:nx);
return;

