% load geomtrace1reg.dat
% load geomtrace2reg.dat
% A=geomtrace1reg;
% B=geomtrace2reg;

load geomtrace1.dat
load geomtrace2.dat
A=geomtrace1;
B=geomtrace2;


sxa=A(:,1);
sya=A(:,2);
gxa=A(:,3);
gya=A(:,4);

sxb=B(:,1);
syb=B(:,2);
gxb=B(:,3);
gyb=B(:,4);

figure(1)
sx=[sxa;sxb];
sy=[sya;syb];
gx=[gxa;gxb];
gy=[gya;gyb];

subplot(221);plot(sx,sy,'o',gx,gy,'+');

title('Acquisition coordinates (2 CDPS)');
legend('shots','receivers');
%vaxis =[-2000 2000 -2000 2000];
%axis(vaxis);

oxa=sxa-gxa;
oya=sya-gya;
oa=sqrt(oxa.^2+oya.^2);

oxb=sxb-gxb;
oyb=syb-gyb;
ob=sqrt(oxb.^2+oyb.^2);

ind=find(oxa==0);oxa(ind)=1;
azimutha=atan(oya./oxa).*180/pi;

ind=find(oxb==0);oxb(ind)=1;
azimuthb=atan(oyb./oxb).*180/pi;

figure(1)
subplot(224);
plot(oa-ob,azimutha-azimuthb,'o');
title('offset vs azimuth (difference)')
%legend('\Delta offset vs \Delta azimuth');
%legend('cdp1','cdp2')


figure(1);
subplot(223);
plot(oxa-oxb,oya-oyb,'o');
title('offsetx vs offsety (difference)')
%legend('offsetx vs offset y');
%legend('cdp1','cdp2')


mxa=(sxa+gxa)/2.;
mya=(sya+gya)/2.;

mxb=(sxb+gxb)/2.;
myb=(syb+gyb)/2.;

figure(1);
subplot(222)
plot(mxa,mya,'o',mxb,myb,'+');
title('midpoint x vs midpoint y');
legend('cdp1','cdp2')
prepfig

