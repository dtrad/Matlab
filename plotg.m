clear
%close all
load reccoord.txt

geom=reccoord;
A=geom(:,:);
scale=1;
sxa=A(:,1)/scale;
sya=A(:,2)/scale;
gxa=A(:,3)/scale;
gya=A(:,4)/scale;
mxa=(sxa+gxa)/2.;
mya=(sya+gya)/2.;

oxa=sxa-gxa;
oya=sya-gya;
oa=sqrt(oxa.^2+oya.^2);


figure(1)
subplot(121);plot(sxa,sya,'.',gxa,gya,'o',mxa,mya,'+');
title('input geom');
legend('shots','receivers','midpoints');


figure(1);
subplot(122);
plot(oxa,oya,'.');
title('offsetx vs offset y');

%print('~/ps/int5dgeombrooks.eps','-deps','-color');


