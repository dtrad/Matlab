function [xx]=shot_gathers_freq2xt(XX)
% [xx]=shot_gathers_freq2xt(XX)
% Given ns shot gathers computes the fft of every shot gather
[nt,ng,ns]=size(XX);

for ii=1:ns
   xx(:,:,ii)=ifft(duplic(XX(:,:,ii)));
end


