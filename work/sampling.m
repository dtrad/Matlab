close all;clear;

NP=512;
FREQ=0.7;
scale=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
te=0:NP-1;
[to,dt]=randdt(NP,0.5);

tso=zeros(1,NP);
tso(round(to)+1)=1;

tse=zeros(1,NP);
tse(round(te(1:10:length(te))+1))=1;

NYQUIST=1/2/max(dt)
NPO=length(tso);
faxiso=(-NPO/2:NPO/2-1)*1/NPO;
faxise=(-NP/2:NP/2-1)*1/NP;

figure,
subplot(221),plot(tse);title('regular');
subplot(222),plot(faxise,abs(fftshift(fft(tse))));
title('fft(reg. sampling function)');

subplot(223),plot(tso);title('irregular');
subplot(224),plot(faxiso,abs(fftshift(fft(tso))));
title('fft(irreg. sampling function)');


x=sin(2*pi*FREQ*te);
xs=sin(2*pi*FREQ*to);

%X=abs(fft(x));
[X,we]=dft(x,te,scale);
[XS,wo]=dft(xs,to,scale);
X=abs(X);
XS=abs(XS);

figure,
subplot(221),plot(te,x);title('Harmonic x, reg sampling'); 
subplot(222),plot(to,xs);title('Harmonic x, irreg sampling');
subplot(223),plot(we/(2*pi),(X));title('dft of x, reg sampling'); 
subplot(224),plot(wo/(2*pi),(XS));title('dft of x, irreg sampling'); 
