function [w]=mywindow(NF,NP,MM)
if nargin<3 MM=4;end
lim1=NF/MM/2;
lim2=NF-NF/MM;
lim3=NF/MM;

wtemp=hanning(lim3);
w=[wtemp(1:lim1);ones(lim2,1);wtemp(lim1+1:lim3)];
w=w*ones(1,NP);
