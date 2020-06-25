%Lab 04
% Amplitude modulation, leakage, spectral resolution, windowing and
% something else.

% Parameters
dt=0.004;
nt=128;
f1=50;
ii=1; % Counter for plots
% Time axis
t=0:nt-1;
t=t*dt;
% Frequency axis
f=freqaxis(dt,nt);
% Signal
x=sin(2*pi*f1*t);

% Long signal
X=fftshift(fft(x));
figure(ii);ii=ii+1;
plotspectra(x,X,t,f,'signal');


% Box car window
z=zeros(size(t));
z(80:100)=1;
Z=fftshift(fft(z));
figure(ii);ii=ii+1;
plotspectra(z,Z,t,f,'box car window');

% Truncated signal or modulated window (frequency shifting)
y=x.*z;
Y=fftshift(fft(y));
figure(ii);ii=ii+1;
plotspectra(y,Y,t,f,'Truncated signal and modulated window');


% Time shifting
t0=0.1;
y2=x.*z;
Y2=fftshift(fft(y2));
Y2=Y2.*exp(i*2*pi*f*t0);
y2=ifft(ifftshift(Y2));
figure(ii);ii=ii+1;
plotspectra(y2,Y2,t,f,'time shifting');

% Hanning window
figure(ii);ii=ii+1;
han=zeros(size(t));
han(80:100)=hanning(21);
HAN=fftshift(fft(han));
plotspectra(han,HAN,t,f,'Hanning window');

% Truncated windowed signal 
y=x.*han;
Y=fftshift(fft(y));
figure(ii);ii=ii+1;
plotspectra(y,Y,t,f,'Hanning.*signal');


