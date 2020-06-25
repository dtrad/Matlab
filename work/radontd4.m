close all
clear
synthrad;
method=1 % Thorson and Claerbout Cm=x^2;
%method=2 % Semblance weigthing Cm=1/(1-s);

plotfig=1; % Plot figure for presentation

iterend=3;
eps1=1e-3;
eps2=1e-3;
inneriter=5;
step=.95;
perc=0.1;
theta=0;
noise=randn(size(d));

[nt np]=size(v);
[nt nh]=size(d);
WW=ones(size(nt*np,1));
dr=radonopd(v,t,h,p,40);
dr=dr+perc/100*noise;
figure,wigb(dr,1,h,t);
if (method==1) 
   WW=ones(nt*np,1);
elseif (method==2)   
   ss=semblance(dr,t,h,p);
   WW=1.1-ss(:);
end  

[vr]=wtcglstdd(t,h,p,WW,dr(:),2,0,step,theta);
figure,wigb(reshape(vr,nt,np));  

  
for i=1:iterend;
   if (method==1) WW=sqrt(1./(eps2+abs(vr(:)))+eps1);
   elseif (method==2) WW=1.1-ss(:);
   end      
	[vr]=wtcglstdd(t,h,p,WW,dr(:),inneriter,0,step,theta);  
   figure,wigb(reshape(vr,nt,np));
   drr=radonopd(vr,t,h,p,theta);
   if (method==2) ss=semblance(drr,t,h,p);end
end
if (plotfig==1)
 drr=radonopd(vr,t,h,p,theta);
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
