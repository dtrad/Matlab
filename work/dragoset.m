% Multiple prediction  (Dragoset - 1982)
% Daniel Trad - UBC
%clear
%close all
[NH,NH,NF]=size(P11);
clear P1;
clear P;V=vel(1);
NSG=NH/2;
[p]=dsw2xt(P11(:,NSG,:));
%p=p.*((t(:)).^0.5*ones(1,NH));
S=S(:);
S=[S(1:100);zeros(2*NF-200,1);S(2*NF-100+1:2*NF)];
SS=S(:)*ones(1,NH);
p=ifft(fft(p).*SS);
p=normalize(p)*coeff(1);
figure,
subplot(121);wigb(p(1:2*NF,:),1,hh-hh(NSG),t);
title('(a) Source gather- Primaries');ylabel('offset');xlabel('time');
k=frequency(dh,NH);
k=k(:);%*ones(1,NH);
for ii=1:NF;
   ww=2*pi*ii/dt/2/NF;
   PRKW=ifft((real(1-(k*V/ww).^2).^0.5).*fft(P11(:,NH/2,ii)));
   P1(:,:,ii)=(1-i)*(ww/4/pi).^0.5*dh*P11(:,:,ii)*PRKW(:);
end;
[p]=dsw2xt(P1(:,1,:));p=-p.*((t(:)*ones(1,NH)).^0.5);
p=ifft(fft(p).*SS);
p=normalize(p)*coeff(1)^2;
subplot(122),wigb(p(1:2*NF,:),1,hh-hh(NSG),t);
title('(b) Source gather- First order multiples');ylabel('offset');xlabel('time');


%p=ifft(fft(p).*SS);
