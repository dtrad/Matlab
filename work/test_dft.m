% Test DFT as a matrix vector multiplication and compare with FFT
t=0:7; % time axis
x=1:8; % time series
X=fft(x)
X2=dft(x,t).'
disp('Maximum difference between DFT and FFT:');
disp(max(abs(X - X2)));

function [X,f]=dfttest(x,t);
% [X,f]=dft(x,t);
% f axis designed to match the fft result exactly
nt=length(x)
dt=(t(end)-t(1))/(length(t)-1);
f=(-nt)/2:(nt-2)/2
f=fftshift(f)
f=f/(dt*nt);
w=f*2*pi;
F=exp(-1i*(w(:)*t(:).'));
X=F*x(:);
return
end