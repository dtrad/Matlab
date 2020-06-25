%close all;
clear;
alfa=1e-5 % hyperparmeter
NP=512;
FREQ1=0.25;
FREQ2=0.7;
scale=1;
A1=1;
A2=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
te=0:NP-1;
[to,dt]=randdt(NP,1);

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

subplot(221),plot(tse);title('(a)');
subplot(222),plot(faxise,abs(fftshift(fft(tse))));title('(b)');
subplot(223),plot(tso);title('(c)');
subplot(224),plot(faxiso,abs(fftshift(fft(tso))));title('(d)');

x=A1*sin(2*pi*FREQ1*te)+A2*sin(2*pi*FREQ2*te);
xs=A1*sin(2*pi*FREQ1*to)+A2*sin(2*pi*FREQ2*to);

%X=abs(fft(x));
[X,we]=dft(x,te,scale,alfa);
[XS,wo]=dft(xs,to,scale,alfa);
%XS=X;wo=we;
%X=XS;we=wo;
X=abs(X);
XS=abs(XS);

figure,
%subplot(221),plot(te,x);title('Harmonic x, reg sampling'); 
%subplot(222),plot(to,xs);title('Harmonic x, irreg sampling');
%subplot(223),plot(we/(2*pi),(X));title('dft of x, reg sampling'); 
%subplot(224),plot(wo/(2*pi),(XS));title('dft of x, irreg sampling'); 

subplot(221),plot(te,x);title('(a)'); 
subplot(222),plot(to,xs);title('(b)');
subplot(223),plot(we/(2*pi),(X));title('(c)'); 
subplot(224),plot(wo/(2*pi),(XS));title('(d)'); 
































