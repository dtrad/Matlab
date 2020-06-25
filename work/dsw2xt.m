function [p]=dsw2xt(P,NF)
if nargin < 2 NF=256;end
P=permute(P,[3 1 2]);
P=zeropadm(P,NF);
%p=ifft(duplic(windowing(P(1:NF,:))));
p=ifft(duplic((P(1:NF,:))));
%figure,wigb(p(1:NF,:),1,hh,tt);
return;


