function [c]=deconv_freq_domain(a,b,lambda)
% [c]=deconv_freq_domain(a,b,lambda)
% Daniel Trad
if nargin<3 lambda=eps;end
la=length(a);
A=fft(a);
b=padzeros(b,la);
B=fft(b);B=B(:);A=A(:);
C=A.*conj(B)./(abs(B).^2+lambda);
c=real(ifft(C));