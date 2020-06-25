function [x,cc,a]=trace(r,V,h,dt,NF)
% x=trace(r,V,h,dt)
% Given the model parameters obtain the trace
% Daniel Trad- UBC - 21-08-98

if nargin < 5 NF=512;end
if nargin < 4 dt=0.004;end
NL=length(V)-1;
c(1)=-1;

tt=0;
for ii=1:NL;
c(ii)=coeff(r(ii),r(ii+1),V(ii),V(ii+1));
dti(ii)=h(ii)/V(ii);
ddt(ii)=t2index(dti(ii),dt);
tt=tt+ddt(ii);
tti(ii)=tt;
end

cc=zeros(NF,1);

cc(1)=c(1);
for ii=1:NL
   cc(tti(ii))=c(ii);
end

% coeff2waves
[x,a]=rf2waves(cc,NF-1);
x(NF)=0;
