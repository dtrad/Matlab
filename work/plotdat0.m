function plotdat0(u,U,ur,axis1)
global dh nh h_near tor t ff  HH dt

figure,
subplot(221),
plotrace(u,dh,nh,h_near,tor);
title('original data in t-x domain'),xlabel('time'),ylabel('offset (meters)'),
axis(axis1);

subplot(222),
temp=abs(U)./max(max(abs(U)));
plotrace(abs(temp),dh,nh,h_near,ff);
title('Data in freq-x domain (spectrum)'),xlabel('freq'),ylabel('offset'),


subplot(223),
plotrace(ur,dh,nh,h_near,t);
title('recovered (freq-time) data in t-x domain'),xlabel('time'),ylabel('offset'),
axis(axis1);

nto=max(size(u));
ntr=max(size(ur));
ntf=min(nto,ntr);

%urtemp=shrinkt2(ur,HH);
res=u(1:ntf,:)-ur(1:ntf,:);%Residuals: 
tt=(0:ntf-1)*0.004;

subplot(224),
plotrace(res,dh,nh,h_near,t);
title('residuals (freq-time) data in t-x domain'),xlabel('time'),ylabel('offset'),
axis(axis1);