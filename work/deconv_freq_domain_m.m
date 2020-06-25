function [c]=deconv_freq_domain_m(a,b,lambda)
% [c]=deconv_freq_domain(a,b)
% Daniel Trad
a=seis_shape(a);
[NTa,NHa]=size(a);
[NTb]=length(b);
%if NTa~=NTb display('number of rows must be equal');break,end
for ii=1:NHa,c(:,ii)=deconv_freq_domain(a(:,ii),b,lambda);end
   
