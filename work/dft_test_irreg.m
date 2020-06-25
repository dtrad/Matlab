function dft_test
dx=50;
kx1=1/(2*dx)*0.05;

xmin=0;
xmax=20000;

% First Test: regular data
xreg=xmin:dx:xmax;nx=length(xreg);
datareg=sin(kx1*2*pi*xreg);
%keyboard;

[DR,kxreg]=dft(datareg,xreg);
DR=fftshift(DR);
DR2=fftshift(fft(datareg));

SPDR=abs(conj(DR).*DR);
max1=max(max(SPDR));

SPDR2=abs(conj(DR2).*DR2);
max2=max(max(SPDR2));


figure,
subplot(211);plot(fftshift(kxreg),SPDR/max1);
title('Regular data: Top DFT, Bottom FFT');
vv=axis;

subplot(212);plot(fftshift(kxreg),SPDR2/max2);
axis(vv);
print -dpng regularDFT.png



% Second Test: irregular data
xi(1)=xmin;
for ix=2:nx
    xi(ix)=xi(ix-1)+dx*rand*2;
end
datairreg=sin(kx1*2*pi*xi);

[D,kxi]=dft(datairreg,xi);
D=fftshift(D);
D2=fftshift(fft(datairreg));

SPD=abs(conj(D).*D);
max1=max(max(SPD));

SPD2=abs(conj(D2).*D2);
max2=max(max(SPD2));

figure;
subplot(211);plot(fftshift(kxi),SPD/max1);
title('Irregular Data: Top DFT, Bottom FFT');
vv=axis;

subplot(212);plot(fftshift(kxi),SPD/max2);
axis(vv);
print -dpng irregularDFT.png

% Probe that the wrong spectrum still recovers data
figure
subplot(211),plot(datairreg);
title('irregular data can be recovered from the wrong spectrum');
subplot(212),plot(ifft(fft(datairreg)));
print -dpng recoverfromwrong.png

figure;
subplot(311);plot(xreg,datareg,'-+');
title('regular data');
vv=axis;

subplot(312);plot(xi,datairreg,'-+');
title('irregular data');
axis(vv);

subplot(313);plot(1:nx,xreg,'o',1:nx,xi,'+');
title('reg vs irreg axis');
print -dpng irregdata.png

% figure;
% subplot(211);plot(kx,atan(imag(D)./real(D)));
% subplot(212);plot(kx,atan(imag(D2)./real(D2)));



function [X,f]=dft(x,t,nz);
% [X,f]=dft(x,t,nz);

if (nargin < 3) nz = 0; end

nt=length(x);
dt=(t(end)-t(1))/(length(t)-1);

ntz=nt+nz;

% original 0 centered axis.
% f=(-ntz+1)/2:(ntz-1)/2;

% generate a frequency axis that matches fft axis.
f=(-ntz+2)/2:(ntz)/2;
vv=fftshift(f);f=[vv(end) vv(1:end-1)];

f=f/(dt*ntz);
w=f*2*pi;

%w=-lw/2:lw/2-1;
%w=w/(lw-1)*2*pi;

F=exp(-i*(w(:)*t(:).'));

X=F*x(:);

return
