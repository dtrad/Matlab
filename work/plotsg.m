[NH,NH,NF]=size(P11);

clear P1;
clear P;
dt=0.004;
dh=15;
hh=0:NH-1;hh=hh*dh;

f=70;
t=0:NF-1;t=t*dt;

wav=rickerm(f,dt);
wav=padzeros(wav,512);
S=fft(wav);
SS=S(:)*ones(1,NH);
figure
for NSG=1:8:32
xa=hh(NSG);
[p]=dsw2xt(P11(:,NSG,:));
%p=p.*((t(:)).^0.5*ones(1,NH));
p=ifft(fft(p).*SS);
subplot(221+round(NSG/8));wigb(p(1:NF,:),1,hh-xa,t);
title('source gathers- Primaries');ylabel('offset');xlabel('time');
end