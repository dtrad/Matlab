function [w]=hamming(n)
w=zeros(n,1);
k=0:n-1;
w(k+1)=0.54-0.46*cos(2*pi*(k/(n-1)));
return