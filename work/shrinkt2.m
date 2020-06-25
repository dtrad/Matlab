function utemp=shrinkt2(u,HH)
% Shrinks a matrix u to a matrix utemp according to vector HH
% Input
% 		u original data
% 		HH vector with zeros and ones: 0 -> skip trace, 1 -> include trace
%
% output 
%   	utemp shrunk data
%
% Daniel Trad-- UBC. Canada. 2-07-98

nh=length(HH);
kk=0;

for ii=1:nh
   if HH(ii)~=0 
      kk=kk+1;
      utemp(:,kk)=u(:,ii);
   end
end;   
   
