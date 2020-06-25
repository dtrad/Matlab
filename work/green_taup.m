
option_xt='y';

GREEN_S=R1.*exp(-i*w*q*(2*z0-z-hs)).*(exp(-2*i*hs*(w*q))-1).*(exp(-2*i*z*w*q)-1);
if option_direct_wave=='y'
   GREEN_INC=exp(-i*w*q*abs(z-hs))-exp(-i*w*q*(z+hs));
   GREEN=i./(2*w*q+eps).*(GREEN_INC+GREEN_S);
else
   GREEN=i./(2*w*q+eps).*(GREEN_S);
end

P0=(Aw(:)*ones(1,NP)).*GREEN;
P0t=ifft(duplic(P0(1:NF/2,:)));

figure(16),
subplot(221),wigb(P0t(1:NF/2,:),1,p(:),ttaxis(1:NF/2));
title('P0(tau,p)+ ghosts'),ylabel('tau(sec)');xlabel('p(s/m)');

GREEN_S=R1.*exp(-i*w*q*(2*z0-z-hs)).*(1-R1.*0.7.*exp(-2*i*z0*w*q)).*(exp(-2*i*hs*(w*q))-1).*(exp(-2*i*z*w*q)-1);
if option_direct_wave=='y'
   GREEN_INC=exp(-i*w*q*abs(z-hs))-exp(-i*w*q*(z+hs));
   GREEN=i./(2*w*q+eps).*(GREEN_INC+GREEN_S);
else
   GREEN=i./(2*w*q+eps).*(GREEN_S);
end

P0=(Aw(:)*ones(1,NP)).*GREEN;
P0t=ifft(duplic(P0(1:NF/2,:)));

figure(16),
subplot(222),wigb(P0t(1:NF/2,:),1,p(:),ttaxis(1:NF/2));
title('P0(tau,p) + ghosts + 1st order mult.')
ylabel('tau(sec)');
xlabel('p(s/m)');

GREEN_S=R1.*exp(-i*w*q*(2*z0-z-hs)).*(1-R1.*0.7.*exp(-2*i*z0*w*q)+(R1.*0.7.*exp(-2*i*z0*w*q)).^2).*...
   (exp(-2*i*hs*(w*q))-1).*(exp(-2*i*z*w*q)-1);
if option_direct_wave=='y'
   GREEN_INC=exp(-i*w*q*abs(z-hs))-exp(-i*w*q*(z+hs));
   GREEN=i./(2*w*q+eps).*(GREEN_INC+GREEN_S);
else
   GREEN=i./(2*w*q+eps).*(GREEN_S);
end

P0=(Aw(:)*ones(1,NP)).*GREEN;
P0t=ifft(duplic(P0(1:NF/2,:)));

figure(16),
subplot(223),wigb(P0t(1:NF/2,:),1,p(:),ttaxis(1:NF/2));
title('P0(tau,p) + ghosts + 1st-2nd order mult.')
ylabel('tau(sec)');
xlabel('p(s/m)');


if option_xt=='y'
p=(0:NP-1)*dp;
h=(0:NH-1)*dh;
P0xt=taup2xt_inv(P0t(1:NF/2,1:NP),10,0.1,'linear',h,p(1:NP),dt);
figure(16),subplot(224),wigb(real(P0xt),1,h,ttaxis(1:NF/2));
ylabel('time (sec)'),xlabel('offset (m)');
mytext=sprintf('Direct Model, Pressure field P%d(t,x) %d layers',ll-1,nl-1)
title(mytext);
end