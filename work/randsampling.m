function [t,dt]=randsampling(lt,dtb,scale)
%     [t,dt]=randsampling(lt,dtb,scale)
% Creates an irregularly sampled axis.
%
% scale controls the degree of irregularity
% dtb, sampling interval base
% lt, length of the output
%
% Daniel Trad

if nargin<3 scale=1.1;end
if nargin<2 dtb=50;end
if nargin<1 lt=63;end

pert=(rand(1,lt))-0.5;
size(pert)

t(1)=0;
dt(1)=0;
for ii=2:lt
   dt(ii)=dtb*(1+scale*pert(ii));
   t(ii)=t(ii-1)+dt(ii);
end



   
