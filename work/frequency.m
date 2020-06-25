function [w]=frequency(dt,NF);
% [w]=frequency(dt,NF);
% Given dt and Number of frequencies returns the angular frequency axis;
% Daniel Trad - UBC
FS=1/dt;
m=0:NF/2-1;
m=m.*FS/NF;
w=2*pi*m';
w(NF/2+1)=pi*FS;
w(NF/2+2:NF)=-w(NF/2:-1:2);
