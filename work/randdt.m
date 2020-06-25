function [t,dt]=randdt(lx,scale)
if nargin<2 scale=0.1;end
pert=(rand(1,lx))-0.5;
dtb=1;
size(pert)

t(1)=0;
dt(1)=0;
for ii=2:lx
   dt(ii)=dtb+scale*pert(ii);
   t(ii)=t(ii-1)+dt(ii);
end


   
