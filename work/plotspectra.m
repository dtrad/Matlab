function plotspectra(x,X,t,f,s);
if (nargin<5) s=[''];end
subplot(311);plot(t,x);title([s,': time domain'])
subplot(312);plot(f,abs(X));title('Amplitude spectrum')
subplot(313);plot(f,unwrap(angle(X)));title('Phase spectrum')

