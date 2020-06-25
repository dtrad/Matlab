function [h1,gap]=h1_without_gaps(h0,dh)
% [gap]=lookfor_gaps(h0,dh)
% Daniel trad_ UBC
h0=h0(:).';
nh=length(h0);
ii=1:nh-1;
dh0=h0(ii)-h0(ii+1);
gap=find(abs(dh0)>2*dh);
ngap=length(gap);
if ngap==0 
   h1=h0;
else   
h1=[h0(1:gap(1))];
for ii=1:ngap-1;
   temp=h0(gap(ii))+dh:dh:h0(gap(ii)+1)-dh/dh;
   temp=temp(:).';
   h1=[h1,temp,h0(gap(ii)+1:gap(ii+1))];
end
temp=h0(gap(ngap))+dh:dh:h0(gap(ngap)+1)-dh/dh;
h1=[h1,temp,h0(gap(ngap)+1:nh)];
end   


