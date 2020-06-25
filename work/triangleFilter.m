function [z]=triangleFilter(nb,nx,x)
% implementation of boxcar twice using recursive filter.
% 
b=zeros(1,nx+nb);
%b(1:nx)=integrate2(x);
b(1:nx)=integrate2fft(x);
%b(1:nx)=cumsum(x(1:nx));
nh=(nb-1)/2;
y=zeros(1,nx+nb);
for (i=nb+1:nx)
 y(i) = 2*b(i) - b(i+nh+1) - b(i-nh-1);
end
scale = 1/((nh+1)*(nh+1));
y=y*scale;
y=y-mean(y);

z(1:nx)=y(1:nx);
return;

