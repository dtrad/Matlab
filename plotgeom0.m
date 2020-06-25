%load wilfoutcdp6006.dat
%B=wilfoutcdp6006;
load datageom.dat

A=datageom;
sxa=A(:,1);
sya=A(:,2);
gxa=A(:,3);
gya=A(:,4);
figure(1)
subplot(221);plot(sxa,sya,'.',gxa,gya,'+');
title('input geom');
legend('shots','receivers');
vaxis =[-2000 2000 -2000 2000];
axis(vaxis);

%select=find(A(:,5)==-135);
%B=A(select,:);

B=A;
sxb=B(:,1);
syb=B(:,2);
gxb=B(:,3);
gyb=B(:,4);

figure(1)
subplot(222);
plot(sxb,syb,'.',gxb,gyb,'+');
title('Output geom');
legend('shots','receivers');
axis(vaxis);

oxa=sxa-gxa;
oya=sya-gya;
oa=sqrt(oxa.^2+oya.^2);

oxb=sxb-gxb;
oyb=syb-gyb;
oa=sqrt(oxb.^2+oyb.^2);

figure(1);
subplot(223);
plot(oxa,oya,'.',oxb,oyb,'+');
title('offsetx vs offset y');
legend('input','output');


mxa=(sxa+gxa)/2.;
mya=(sya+gya)/2.;

mxb=(sxb+gxb)/2.;
myb=(syb+gyb)/2.;

figure(1);
subplot(224)
plot(mxa,mya,'.',mxb,myb,'+');
title('midpoint x vs midpoint y');
legend('input','output');

%print('~/ps/int5dgeomcomp.eps','-deps','-color');


