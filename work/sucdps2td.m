close all
clear
load cdps2.mat 
method=1 % Thorson and Claerbout Cm=x^2;
method=2 % Semblance weigthing Cm=1/(1-s);
plotfig=1; % Plot figure for presentation

iterend=3;
eps1=1e-3;
eps2=1e-3;
inneriter=5;
step=.9;
perc=20;

noise=rand(size(d))-.5;

np=100;
pmin=0;
pmax=3.4e-7;
dp=(pmax-pmin)/(np-1);
p=pmin:dp:pmax;

[nt nh]=size(d);


figure,wigb(d,1,h,t);
if (method==1) 
   WW=ones(nt*np,1);
elseif (method==2)   
   ss=semblance(d,t,h,p);
   WW=1.1-ss(:);
end  

[vr]=wtcglsm(d(:),t,h,p,WW,2,0,step);
figure,wigb(reshape(vr,nt,np));  
  
for i=1:iterend;
   if (method==1) WW=(1./(eps2+abs(vr(:)))+eps1);
   elseif (method==2) WW=1.1-ss(:);
   end      
   [vr]=wtcglsm(d(:),t,h,p,WW,inneriter,0,step);  
   figure,wigb(reshape(vr,nt,np));
   drr=radonop(vr,t,h,p,WW);
   if (method==2) ss=semblance(drr,t,h,p);end
end
if (plotfig==1)
 drr=radonop(vr,t,h,p,WW);
 figure,
 subplot(221);wigb(reshape(vr,nt,np),1,p,t);
 %subplot(223);wigb(reshape(vr,nt,np),1,p,t);
 subplot(222);wigb(reshape(d,nt,nh),1,h,t);
 subplot(223);wigb(reshape(drr,nt,nh),1,h,t);

 subplot(221),ylabel('time (s)');
 subplot(222),ylabel('time (s)');
 %subplot(223),ylabel('time (s)');
 subplot(223),ylabel('time (s)');
 subplot(221),title('q (s/m)^2')
 subplot(222),title('offset (m)')
 subplot(221);xlabel('(a)')
 subplot(222);xlabel('(b)')
 %subplot(223);xlabel('(c)')
 subplot(223);xlabel('(c)')
end






