function [p]=hyperb_scales(p,sx,h,st,nt,nx,dt,dx)
p=zeros(nt,nx);
p=0
for ix=nx/2:nx
  x=(ix-nx/2)*dx,
  t=0:nt-1;t*dt;
  for i=1:nt
    tmp=((t(i)-sqrt((x/sx)^2+h*h))/st)^2;
    p(i,ix)=((1-tmp)*exp(-tmp/2))
  end
end
    