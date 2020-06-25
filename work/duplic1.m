function [out]=duplic1(in)
% Given half of the FFT series duplic.m  duplicates the data to recover 
% the full spectrum such that the original time series is real.
% the input can be a single vector or a matrix, with the frequency dimension
% longer than the spatial dimension
%
% [out]=duplic(in)
%
%  Daniel Trad, UBC- 30/03/98

%[NF NP]=size(in);if NF<NP in=in.';end
[NF NP]=size(in);
 m=1:NF;
 j=1:NP;
 out(m,j)=in(m,j);

% Nyquist Freq is not in the data so that it must be set 
% either = 0 or = previous value.
 out(NF+1,j)=zeros(size(j));  % Nyquist

 %out(NF+1,j)=real(out(NF,j));
 %out(NF+1,j)=(out(NF,j));   
 
 m=1:NF-1;
 out(NF+1+m,j)=conj(in(NF+1-m,j));
 
 %size(out)
 %size(in)

