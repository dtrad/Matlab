function [S,freq] = stransform(x,t,fmin,fmax,n);
%
% STRANSFORM Syntax:  [S,f] = stransform(x,t,fmin,fmax,n)
%
% where:  S = output time frequency representation
%         f = output frequency vector (MHz)
%         x = signal to be analyzed
%         t = time vector (ns)
%         fmin = minimum frequency to analyze in MHz (must be positive)
%         fmax = maximum frequency to analyze in MHz (must be <= Nyquist)
%         n = number of points to use for fft
%
% This function computes the S-transform time-frequency representation of a signal,
% which is a hybrid STFT/WT, using the fast frequency domain algorithm provided in 
% Stockwell et al. (1996).
%
% by James Irving
% June, 2001

t = t(:).';        % convert t and x to row vectors if necessary
x = x(:).';
si = t(2)-t(1);    % sampling interval (s)
nppt = length(t); % number of points per trace
fn = 1/(2*si);     % Nyquist frequency (Hz)

if (nargin < 5)    n = 2^nextpow2(nppt); end

df = 2*fn/n;       % sampling interval in frequency domain (Hz)

if (nargin < 4)    fmin=0; end

if (nargin < 3);  fmax=fn; end

 
f = -fn:df:fn-df;  % frequency vector (Hz) for FFT stuff
w = 2*pi*f;        % frequency vector (rad/s) for FFT stuff
f2 = 0:df:fn;      % frequency vector (Hz) for computing S-transform frequencies
w2 = 2*pi*f2;      % frequency vector (rad/s) for computing S-transform frequencies

fmin = fmin/1e3;                   % minimum analysis frequency (Hz)
fmax = fmax/1e3;                   % maximum analysis frequency (Hz)
[temp,jmin] = min(abs(f2-fmin));   % closest index of min. frequency in f2
[temp,jmax] = min(abs(f2-fmax));   % closest index of max. frequency in f2

X = fftshift(fft(x,n));            % Fourier transform input trace
X = [X,X];                         % double to handle wrap around

for j=jmin:jmax                    % perform S-transform
   if j==1 
      S(1,:) = mean(x)*ones(size(x));  % set DC value of transform to average of signal
   else 
      G = exp(-2*pi^2*w2(j)^-2.*w.^2); % Gaussian window
      H = X(j:j+n-1);                  % shifted Fourier spectrum
      B = H.*G;                        % fft of S-transform for frequency f2(j)
      temp = ifft(fftshift(B));        % tranform to time
   	S(j-jmin+1,:) = temp(1:nppt);    % place in matrix S
   end
end

freq = f2(jmin:jmax)*1e3;          % output frequency vector (MHz)
