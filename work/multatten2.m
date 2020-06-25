% Multiple prediction  (Dragoset - 1982)
% Daniel Trad - UBC
clear
close all

%load c:\daniel\thesis\P00;P00=P11;clear P11;
load P11;
[NH,NH,NF]=size(P11);
R1=1*eye(NH);
NSG=32;
NNUL=NF-100;
dt=0.004;t=(0:2*NF-1)*dt;
dh=15;hh=0:NH-1;hh=hh*dh;
f=70;wav=ricker(f,dt);wav=padzeros(wav,2*NF);
S=fft(wav);
coeff(1)=0.2

[p1]=dsw2xt(P11(:,NSG,:));

S=S(:);
S=[S(1:100);zeros(2*NF-200,1);S(2*NF-100+1:2*NF)];
SS=S(:)*ones(1,NH);
%AW=abs(permute(P11(1,1,:),[3 1 2]));
p1=ifft(fft(p1).*SS);
%p1=normalize(p1)*coeff(1);
figure(1),
subplot(121);wigb(p1(1:2*NF,:),1,hh-hh(NSG),t);
title('Data (a)');ylabel('offset');xlabel('time');
niter=1;
for ii=1:NF-NNUL;
   P1=P11(:,:,ii);
   %P0=P00(:,:,ii);
   P10=P1;
   for jj=1:niter;
      %P10=P1-P10*R1*P1;
      AMP=1./max(max(abs(P1)));
      %AMP=-1./abs(P11(:,:,ii)+1e-3);
      %AMP=-1./2;
      P10=-R1*(AMP*(P10.*P1));
   end
   P11(:,:,ii)=P10;
end;
P11(:,:,NF-NNUL+1:NF)=0;
[p]=dsw2xt(P11(:,NSG,:));
p=ifft(fft(p).*SS);
%p=normalize(p)*coeff(1);
subplot(122),wigb(p(1:2*NF,:),1,hh-hh(NSG),t);
title('(b)');xlabel('offset');ylabel('time');

figure(2),
subplot(121),wigb(p1,1,hh-hh(NSG),t);
title('(a)');xlabel('offset');ylabel('time');
subplot(122),wigb(p1+p,1,hh-hh(NSG),t);
title('(b)');xlabel('offset');ylabel('time');


figure(3),
subplot(121),wigb(p1(:,10:14));
subplot(122),wigb(p(:,10:14));
