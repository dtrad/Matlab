function [z]=mix_fft(x,y)
% [z]=mix_fft(x,y)
% Given two time series x,y  it computes z with amplitude spectrum from y and phase from x.
% Daniel Trad - 14-10-98 UBC
%x=[x(:);zeros(size(x(:)))];
%y=[y(:);zeros(size(y(:)))];
[NF,NH]=size(x);
ampx=abs(fft(x));
phasex=angle(fft(x));
ampy=abs(fft(y));
phasey=angle(fft(y));
z=ifft(ampy.*exp(i*(phasex)));

%z=z(1:length(z)/2);