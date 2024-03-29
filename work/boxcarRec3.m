function [y]=boxcarRec3(nb,nx,x)

% implementation of boxcar using recursive filter.
% replace integral by 1./|omega|
b=zeros(1,nx+nb);
b(1:nx)=cumsum(x(1:nx));
nh=(nb-1)/2;
for (i=nh+1:nx)
    y(i) = b(i+nh) - b(i-nh);
end
y=y/nb*2;

return;