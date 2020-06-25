% Resampling using DFT

close all;clear;

NT=128;
FREQ1=5;
FREQ2=20;
dt=0.004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr=0:NT-1;tr=tr*dt;
[ti,dti]=randsampling(NT,dt,2);

f=freqaxis(dt,NT);


tsi=zeros(1,NT);
tsi(round(ti/dt)+1)=1;
tsi=tsi(1:NT);
tsr=zeros(1,NT);
tsr(round(tr/dt+1))=1;


figure,
subplot(221),plot(tsr);title('regular');
subplot(222),plot(f,abs(fftshift(fft(tsr))));
title('fft(reg. sampling function)');

subplot(223),plot(tsi);title('irregular');
subplot(224),plot(f,abs(fftshift(fft(tsi))));
title('fft(irreg. sampling function)');






xr=sin(2*pi*FREQ1*tr)+cos(2*pi*FREQ2*tr);
xi=sin(2*pi*FREQ1*ti)+cos(2*pi*FREQ2*ti);

xr=xr(:);
xi=xi(:);

figure(1)
subplot(211),plot(tr,xr,'.');v=axis;
subplot(212),plot(ti,xi,'.');axis(v);


figure(2)
subplot(211),plot(f,fftshift(abs(fft(xr))));v=axis;
subplot(212),plot(f,fftshift(abs(fft(xi))));axis(v);

w=2*pi*f;

method='inv'
method='dft'
method='spa'
method='pap'

if (method=='inv')
  XR=dft(xr,tr,w,NT/2);
  XI=dft(xi,ti,w,NT/2);
elseif (method=='dft') %DFT
  XR=dft(xr,tr,w,1);
  XI=dft(xi,ti,w,1);
elseif (method=='spa') %DFT
  XR=dft(xr,tr,w,0);
  XI=dft(xi,ti,w,0);
elseif (method=='pap')
  XR=dft(xr,tr,w,1);
  XIdft=dft(xi,ti,w,1);
  XIinv=dft(xi,ti,w,NT/2);
  XIspa=dft(xi,ti,w,0);
  figure(10)
  v=[min(f) max(f) 0 0.5];
  subplot(411),plot(f,abs(XR));axis(v);title(['(a) FFT from regular' ...
		    ' sampled data'])
  subplot(412),plot(f,abs(XIdft));axis(v);title('(b) DFT');
  subplot(413),plot(f,abs(XIinv));axis(v);title('(c) LSFT');
  subplot(414),plot(f,abs(XIspa));axis(v);title('(d) sparse LSFT');
else
  display ('option not implemented')
end

if (method=='pap')
  xpr=real(idft(XR,tr,w));
  xpirdft=real(idft(XIdft,tr,w));
  xpirinv=real(idft(XIinv,tr,w));
  xpirspa=real(idft(XIspa,tr,w));
  figure(11)
  v=[0 0.5 -2 2];
  subplot(411),plot(tr,xpr,'.',tr,xr);axis(v);title(['(a) Prediction' ...
		    ' from  FFT (reg sampled data)'])
  subplot(412),plot(ti,xpirdft,'.',tr,xr);axis(v);title(['(b)' ...
		    ' Prediction by DFT']);
  subplot(413),plot(tr,xpirinv,'.',tr,xr);axis(v);title(['(c)' ...
		    ' Prediction by LSFT']);
  subplot(414),plot(tr,xpirspa,'.',tr,xr);axis(v);title(['(d)' ...
		    ' Prediction sparse LSFT']);
  
  
else

  XI=XI.*hanning(NT);


  figure(3)
  subplot(211),plot(f,abs(XR));v=axis;
  subplot(212),plot(f,abs(XI));axis(v);

  xpr=real(idft(XR,tr,w));
  xpi=real(idft(XI,ti,w));
  xpir=real(idft(XI,tr,w));


  figure(4)
  subplot(311),plot(tr,xpr,'.');v=axis;title(method)
  subplot(312),plot(ti,xpi,'.');axis(v);
  subplot(313),plot(tr,xpir,'.');axis(v);
  
  figure(5)
  subplot(311),plot(tr,xpr,'.',tr,xr);v=axis;title(method)
  subplot(312),plot(ti,xpi,'.',tr,xr);axis(v);
  subplot(313),plot(tr,xpir,'.',tr,xr);axis(v);
end



