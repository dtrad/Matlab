function plotsvdvectors3(F,W,v,d,data,datart,h,p,t,datartls,datarthr)
% Plot the singular vectors for the kernel and for the
% preconditioned kernel
% Input 
%      F kernel
%      W model weights
% Paper example
%  load data
%  [datartls,datathr]=testradon(data,h,q,t); 
%  load temp  (This file is save by radonf1 when f==25)
%  plotsvdvectors3(F,WW,v,d,data,datart,h,q,t,datartls,datarthr);

[n,m]=size(F);
iid=1:n;
iim=1:m;
[U,S,V]=svd(F);


figure(1),
subplot(421)
wigb(data,1,h,t);title('Offset(m)');xlabel('(a)');ylabel('Time(s)')

subplot(422)
wigb(datart,1,iim,t);title('q(s/m^2)');xlabel('(b)');ylabel('Time(s)')


subplot(423)
plot(iid,real(d));ylabel('Real(m(\omega_i,q))');xlabel('(c)');title('q(s/m^2)')

subplot(425)
plot(iid,real(U(:,3:3)));ylabel('Real(m(\omega_i,q))');xlabel('(d)');title('q(s/m^2)')

FWWI=F*inv(W);
[U,S,V]=svd(FWWI);

subplot(426)
plot(iid,real(U(:,3:3)));;
ylabel('Real(m(\omega_i,q))');xlabel('(e)');title('q(s/m^2)')

subplot(427)
wigb(datartls,1,iim,t);
title('q(s/m^2)');xlabel('(f)');ylabel('Time(s)')

subplot(428)
wigb(datarthr,1,iim,t);
title('q(s/m^2)');xlabel('(g)');ylabel('Time(s)')



