function [S2,freq,xi] = stransform_filter(x,t,dfx,fmin,fmax,n);
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

if nargin==0; help stransform;return;end
   
nt = length(t);  % number of points per trace
dt = t(2)-t(1);    % sampling interval (s)
fn = 1/(2*dt);     % Nyquist frequency (Hz)

if nargin<6; n = 2^nextpow2(nt), end
if nargin<5; fmax=fn, end
if nargin<4; fmin=0, end

t = t(:).';        % convert t and x to row vectors if necessary
x = x(:).';
if (length(x)~=length(t)); 
  display('error: should set equal lengths for x and t');
  return;
end
df = 2*fn/n;       % sampling interval in frequency domain (Hz)
f = -fn:df:fn-df;  % frequency vector (Hz) for FFT stuff
w = 2*pi*f;        % frequency vector (rad/s) for FFT stuff
f2 = 0:df:fn;      % frequency vector (Hz) for computing S-transform frequencies
w2 = 2*pi*f2;      % frequency vector (rad/s) for computing S-transform frequencies

%fmin = fmin/1e3;                   % minimum analysis frequency (Hz)
%fmax = fmax/1e3;                   % maximum analysis frequency (Hz)
[temp,jmin] = min(abs(f2-fmin)),   % closest index of min. frequency in f2
[temp,jmax] = min(abs(f2-fmax)),   % closest index of max. frequency in f2

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
      S(j-jmin+1,:) = temp(1:nt);    % place in matrix S

   end
end

% Shift frequencies up by df. This shift is time and frequency
% dependent. 
%keyboard
S2=S;
fit=50;
lastj=round(jmax/2);
for it=fit:nt
  for j=lastj:-1:jmin
    dfshift=dfx(it)*f2(j);
    idf=round((dfshift)/df);
    if ( (j-jmin+1-idf) > jmin )
      %[idf, j, it]
      S2(j-jmin+1,it)=S(j-jmin+1,it)+0*S(j-jmin+1-idf,it);
    end
  end
end

  
freq = f2(jmin:jmax);          % output frequency vector (MHz)

% Inverse S transform
SS=sum(S2,2);
xi=real(ifft(duplic(SS)));

