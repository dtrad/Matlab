function [rw,t] = ricker_2nd(f,dt,tlength,t0)

%RICKER creates an second order ricker wavelet signal
%
%   RICKER creates and plots a default causal ricker wavelet with:
%
%       peak frequency   = 20 Hz
%       sampling time    = 0.001 seconds
%       length of signal in seconds = 3/f;
%       peak location    = 1/F = 1/20Hz
%
%   RW = RICKER(...) returns the default wavelet in RW.
%
%   [RW,T] = RICKER(...) returns the time vector in T.
%
%   Specifying parameters of peak frequency (F, Hz), number of points (N),
%   and sampling time (DT) are specified by the syntax:
%
%       [RW,T] = RICKER(F)
%       [RW,T] = RICKER(F,DT)
%       [RW,T] = RICKER(F,DT,TLENGTH)
%       
%   [RW,T] = RICKER(F,DT,TLENGTH,T0) creates a ricker wavelet with peak centered
%   at T0.

%   Example :
%    % create a ricker wavelet with 40 Hz 0.02,and 0.1 s length  s between
%    % samples
%    [rw,t] = ricker(40,0.002,0.1);
%    plot(t,rw), xlabel('Time'), ylabel('Amplitude')

% Define inputs if needed
switch nargin
    case 0
        f  = 20;
        dt = 0.001;  
        n = 100;
        t0 = 1/f;
    case 1
        dt = 0.001;  
        n = 100;
        t0 = 1/f;
    case 2
        tlength=3/f;
        n  = fix(tlength/dt)+1;
        t0 = 1/f;
    case 3
        n  = fix(tlength/dt)+1;
        t0 = 1/f;  
   case 4
        n  = fix(tlength/dt)+1;
    otherwise
        warning('RICKER:tooManyInputs','Ignoring additional inputs')
end

% Create the wavelet and shift in time if needed
T = dt*(n-1);
t = 0:dt:T;
tau = t-t0;
a=pi*tau*f;
s = normalize((1-2*a.^2).*exp(-a.^2));

    
 
if nargout == 0
    figure;
        plot(t,s)
        xlabel('Time (s)')
        ylabel('Amplitude')
        title('Ricker Wavelet')
else
        rw = s;
end