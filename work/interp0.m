load /home/dtrad/matlab/xxc.mat;
z=seis_shape(xxc(1:256,1:40));
z=normalize(z);

select=[ones(1,10) ones(1,10) zeros(1,3) ones(1,20-3) zeros(1,10)];
x=shrinkt2(z,select);clear z;
temp=1:length(select);temp=temp*10;
h0=shrinkt2(temp,select);clear temp;

[nt,nh]=size(x);
dh=10;

%temp1=(0:9)*dh;
%temp2=(10:19)*dh; 
%temp3=(20:39)*dh;
%temp4=(40:50)*dh;

%h0=[temp1 temp2 temp3];
%h1=[temp1 temp2 temp3 temp4];

[h1,gap]=h1_without_gaps(h0,dh);
lh1=length(h1);
h1=[h1,h1(lh1)+dh:dh:h1(lh1)+5*dh];

dt=0.004;
pnorm=1;
option='Hampson';
vel_min=500;
np=80;dp=1/vel_min/np;
niter=5;
hyperparameter=1e-1;

[u,haxis,ttaxis,vr,ppaxis]=interp_taup1(x,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter);
figure,wigb(vr,1,ppaxis,ttaxis);
%[vel,pr]=alfa2vel(ppaxis,[0 dp]);
%pr=padzeros(pr,length(ppaxis));
%vel(1:10)=0;
%subplot(211),plot(vel);title('Vel');xlabel('#trace')
%subplot(212),plot(ppaxis);title('alfa');title('alfa');xlabel('#trace')




 
