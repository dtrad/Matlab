p=(0:NP-1)*dp;
h=(0:NH-1)*dh;
P0xt=taup2xt_inv(P0t(1:NF/2,1:40),10,0.1,'linear',h,p(1:40),dt);
figure,subplot(111),wigb(real(P0xt),1,h,ttaxis(1:NF/2));
ylabel('time (sec)'),xlabel('offset (m)');
mytext=sprintf('Direct Model, Pressure field P%d(t,x) %d layers',ll,nl)
title(mytext);
