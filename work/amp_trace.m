function [aaa,amp,dist]=amp_trace(dt,dh,aa,vel,sx,sy)
% [aaa]=amp_trace(dt,dh,aa,vel,sx,sy)
aa=normalize(aa);
[nt,nh]=size(aa);
yy=1:nt;yy=yy(:)*ones(1,nh);
xx=1:nh;xx=ones(nt,1)*xx(:).';
dist=abs(sqrt((vel*dt)^2*(yy-sy).^2+dh^2*(xx-sx).^2));
amp=(1+dist);
aaa=amp.*aa;
