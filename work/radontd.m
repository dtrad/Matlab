close all
clear
synthrad;

plotfig=1; % Plot figure for presentation

iterend=3;
eps1=1e-3;
eps2=1e-3;
inneriter=5;
step=.9;
perc=0;

noise=rand(size(d))-.5;

[nt np]=size(v);
[nt nh]=size(d);

dr=radonop(v,t,h,p);
dr=dr+perc/100*noise;
figure,wigb(dr,1,h,t);
WW=ones(nt*np,1);
[vr]=wtcglstd(t,h,p,WW,dr(:),2,0,1);
figure,wigb(reshape(vr,nt,np));  
  
for i=1:iterend;
WW=(1./(eps2+abs(vr(:)))+eps1);
%WW=(1./(eps2+(vr(:).^2))+eps1);

[vr]=wtcglstd(t,h,p,WW,dr(:),inneriter,0,step);  

figure,wigb(reshape(vr,nt,np));  
end
if (plotfig==1)
 drr=radonop(vr,t,h,p);
 figure,
 subplot(221);wigb(reshape(v,nt,np),1,p,t);
 subplot(223);wigb(reshape(vr,nt,np),1,p,t);
 subplot(222);wigb(reshape(dr,nt,nh),1,h,t);
 subplot(224);wigb(reshape(drr,nt,nh),1,h,t);

 subplot(221),ylabel('time (s)');
 subplot(222),ylabel('time (s)');
 subplot(223),ylabel('time (s)');
 subplot(224),ylabel('time (s)');
 subplot(221),title('q (s/m)^2')
 subplot(222),title('offset (m)')
 subplot(221);xlabel('(a)')
 subplot(222);xlabel('(b)')
 subplot(223);xlabel('(c)')
 subplot(224);xlabel('(d)')
end