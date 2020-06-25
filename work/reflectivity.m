function [cc]=reflectivity(c,lx,threshold)
if nargin<3 threshold=0.00001;end
cc=random('Normal',0,0.13,1,lx);
indc=find(c~=0);
indcc=find(abs(cc)>=threshold);
ccw=zeros(1,lx);
ccw(indcc)=cc(indcc);ccw=padzeros(ccw,lx);
ccw(indc)=c(indc);
cc=ccw;