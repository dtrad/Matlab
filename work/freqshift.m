dt=0.004;
nt=512;
freq=20;
t0=1;
t=[0:nt-1]'*0.004;
    
w=ricker(freq,dt);
x=zeropadm(w,nt);
f=freqaxis(dt,length(x))';
x2=real(ifft(fft(x).*exp(-i*f*2*pi*t0)));
subplot(211);plot(t,x);
subplot(212);plot(t,x2);