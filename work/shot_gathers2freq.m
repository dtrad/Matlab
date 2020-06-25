function [XX]=shot_gathers2freq(xx)
% [XX]=shot_gathers2freq(xx)
% Given ns shot gathers computes the fft of every shot gather
[nt,ng,ns]=size(xx);
%w=hanning(nt)*ones(1,ng);

for ii=1:ns
   XX(:,:,ii)=fft(xx(:,:,ii));
end


