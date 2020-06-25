function [S,f] = stransform2(x,t,num,fmin,fmax);
%
% STRANSFORM Syntax:  [S,f] = stransform(x,t,num,fmin,fmax)
%
% where:  S = output time frequency representation
%         f = output frequency vector (MHz)
%         x = signal to be analyzed
%         t = time vector (ns)
%         num = number of points to use in frequency
%         fmin = minimum frequency to analyze in MHz
%         fmax = maximum frequency to analyze in MHz (should be less than Nyquist)
%
% This function computes the S-transform time-frequency representation of a signal,
% which is a hybrid STFT/WT (see Geophysics v.65 n.4, p.1333).  Watch out for end
% effects.
%
% by James Irving
% October, 2000

if size(t,1)~=1; t=t'; end  % convert t and x to row vectors if necessary
if size(x,1)~=1; x=x'; end

nppt = length(t);           % number of points per trace
%si = t(2)-t(1);             % sampling interval (s)
%fn = 1/(2*si);              % Nyquist frequency (Hz)

if nargin==3; fmin=0; fmax=fn; end
if nargin==4; fmax=fn; end

f = linspace(fmin,fmax,num);   % frequency vector (Hz)
f = f(:);                        % convert to column vector
KF=exp(-i*2*pi*f*t);
const=abs(f)/sqrt(2*pi);

for j=1:nppt
   %disp(j)
   %KT=exp(-0.5*(f*(t(j)-t)).^2);   
   for k=1:num
      %window = const(k)*(KT(k,:)).*KF(k,:);
      window = const(k)*exp(-0.5*(t(j)-t).^2*f(k)^2).*KF(k,:);
      S(k,j) = sum(x.*window);
   end
end

S = S./max(max(abs(S)));


