function [wavl,Aw]=seis_wavelet(freq,NF,dt)
% [wavl,Aw]=seis_wavelet(freq,dt,NF)
if nargin < 3 dt=0.004;end
if nargin < 2 NF=512;end
if nargin < 1 freq=50;end
wav=rickerm(freq,dt);
wavl=padzeros(wav,NF);
Aw=fft(wavl);
