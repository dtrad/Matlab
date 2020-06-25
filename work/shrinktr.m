function [utemp]=shrinktr(u,HH)
% [utemp]=shrinktr(u,HH)
% Select the traces to use. 
% This function is for trying interpolation.
% Input
% 		xxc original data
% output 
%   	u are the data to use
% Daniel Trad-- UBC. Canada. 2-07-98

nh=max(size(HH));
kk=0;
for ii=1:nh
   if HH(ii)~=0 
      kk=kk+1;
      utemp(kk)=u(ii);
   end
end;   
   
