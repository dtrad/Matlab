function [h,p]=radonaxis(dh,nh,np,vmin,vmax)
h=0:nh-1;
h=h*dh;
pmin=1/(vmax^2);
pmax=1/(vmin^2);
dp=(pmax-pmin)/(np-1);
p=pmin:dp:pmax;
