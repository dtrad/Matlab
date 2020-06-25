function [w_min] = min_phase_wavelet(w);  
%
%MIN_PHASE_WAVELET  Given an arbitrary wavelet
%                   this codes retrieves the minimum
%                   phase wavelet using the Hilbert transform
%
%  [w_min] = MIN_PHASE_WAVELET(w)
%
%  IN.
%   w: a wavelet of arbitrary phase (1 column or 1 row)
%
%  OUT.
%   w_min: a min phase wavelet 
%
%  Example 1:
%
%   w = [1,2];     % a maximun phase wavelet
%   w_min = min_phase_wavelet(w);  
% 
%  Example 2:
%
%   w = ricker(0.004,30);    % A ricker of 30Hz sampled every 4msec
%   w_min = min_phase_wavelet(w);
%   figure(1); plot(w); hold on; plot(w_min)
%
%  Example 3:
%   
%   w_min = min_phase_wavelet(ricker(0.004,30)); 
%
% 
%  M.D.Sacchi, July 1998, Dept. of Physics, UofA.
%  sacchi@phys.ualberta.ca


nw = max(size(w));    % lenght of the wavelet
nfft = 4*( 2^nextpow2(nw));

W = log ( abs(fft(w,nfft)) );
W = ifft(W);

for i=nfft/2+2:nfft; W(i)=0.;end;
W = 2.*W;
W(1) =W(1)/2.;

W = exp(fft(W)) ;
w_min = real(ifft(W));
w_min = w_min(1:nw); 

return

