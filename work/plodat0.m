function plotdat0(u,U,ur)
global dh nh h_near tor t ff  

figure,
subplot(221),
plotrace(u,dh,nh,h_near,tor);
title('original data in t-x domain'),xlabel('time'),ylabel('offset (meters)'),

subplot(222),
temp=abs(U)./max(max(abs(U)));
plotrace(abs(temp),dh,nh,h_near,ff);
title('Data in freq-x domain (spectrum)'),xlabel('freq'),ylabel('offset'),

subplot(223),
plotrace(ur,dh,nh,h_near,t);
title('recovered (freq-time) data in t-x domain'),xlabel('time'),ylabel('offset'),

nto=max(size(ur));

res=u-ur(1:nto,:);%Residuals: 
subplot(224),
plotrace(res,dh,nh,h_near,tor);
title('residuals (freq-time) data in t-x domain'),xlabel('time'),ylabel('offset'),
