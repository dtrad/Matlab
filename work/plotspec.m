function plotspec(wavelet,dt,nspec)
% function plotspec(wavelet,dt,nspec)
% Daniel Trad - UBC 2000.
if (nargin < 3) nspec=length(wavelet);end

wavelet=[wavelet(:);zeros(nspec-length(wavelet),1)];
plot(freqaxis(dt,length(wavelet)),fftshift(abs(fft(wavelet))))
title('amplitude spectrum');
xlabel('frequency in Hertz')
ylabel('Amp');