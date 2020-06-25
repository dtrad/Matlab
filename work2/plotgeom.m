load datageomireg.dat
A=datageomireg;
load datageom.dat
A=datageom;
% a=find((A(:,1)<600)&(A(:,1)>400));
% b=find((A(:,2)<600)&(A(:,2)>400));
% c=find((A(:,3)<600)&(A(:,3)>400));
% d=find((A(:,4)<600)&(A(:,4)>400));
% 
% sxa=A(a,1);
% sya=A(b,2);
% gxa=A(c,3);
% gya=A(d,4);

sxa=A(:,1);
sya=A(:,2);
gxa=A(:,3);
gya=A(:,4);
figure(1)
subplot(221);plot(sxa,sya,'.',gxa,gya,'.');
title('Acquisition coordinates');
legend('shots','receivers');
vaxis =[0 1000 0 1000];
axis(vaxis);

oxa=sxa-gxa;
oya=sya-gya;
oa=sqrt(oxa.^2+oya.^2);

ind=find(oxa==0);oxa(ind)=1;
azimuth=atan(oya./oxa).*180./pi;

figure(1)
subplot(224);
plot(oa,azimuth,'.');
title('processing coordinates 2b')
legend('offset vs azimuth');
vaxis =[350 500 -35 35];
axis(vaxis);

figure(1);
subplot(223);
plot(oxa,oya,'.');
title('processing coordinates 2a')
legend('offsetx vs offset y');
vaxis =[-100 100 -100 100];
axis(vaxis);

mxa=(sxa+gxa)/2.;
mya=(sya+gya)/2.;

figure(1);
subplot(222)
plot(mxa,mya,'.');
legend('midpoint x vs midpoint y');
title('processing coordinates 1');
vaxis =[350 450 350 450];
axis(vaxis);

prepfig
%print('~/ps/int5dgeomcomp.eps','-deps','-color');


