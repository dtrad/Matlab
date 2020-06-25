function [S,s] = stransf0(x);
%
% STRANSFORM Syntax:  [S,f] = stransform(x,t,fmin,fmax,n)
%
% where:  S = output time frequency representation
%         x = signal to be analyzed
%
% This function computes the S-transform time-frequency representation of a signal,
% which is a hybrid STFT/WT, using the fast frequency domain algorithm provided in 
% Stockwell et al. (1996).
%
% by James Irving
% June, 2001
% Modified by Daniel Trad 
% October, 2001


N=2^nextpow2(length(x));

if nargin==0; display('usage: [S]=stransf0(x)');end


X = fft(x,N);            % Fourier transform input trace
X = [X,X];                         % double to handle wrap around

m=0:N-1;
n=m;
j=m;
S=zeros(N,N);
for ii=1:N                    % perform S-transform, voice by voice 
   if ii==1 
      S(1,:) = mean(x)*ones(1,N);  % set DC value of transform to average of signal
   else 
      G = exp(-2*pi^2*m.^2./(n(ii))^2); % Gaussian window
      H = X(ii:ii+N-1);     % shifted Fourier spectrum= FT(h(t).exp(-2pift))
      B = H.*G;             % Convolution p(t)*g(t) where p(t)=h(t)exp(-i2pift)
      temp = ifft((B));        % tranform to time
      S(ii,:) = temp(1:N);    % place in matrix S
   end
end

% Inverse ifft;
s=ifft(sum(S,2));