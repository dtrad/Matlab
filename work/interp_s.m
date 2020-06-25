% Load data generated with synth0. xxc contains dh,dt and h_near
load /home/dtrad/matlab/xxc;
h1=h;
z=seis_shape(xxc(1:512,1:101));clear data;
z=normalize(z);
%dh=25;
select=[ones(1,41) zeros(1,20) ones(1,40)];
x=shrinkt2(z,select);clear z;
temp=0:length(select)-1;temp=temp*dh;
h0=shrinkt2(temp,select);clear temp;

[nt,nh]=size(x);

%temp1=(0:9)*dh;
%temp2=(10:19)*dh; 
%temp3=(20:39)*dh;
%temp4=(40:50)*dh;

%h0=[temp1 temp2 temp3];
%h1=[temp1 temp2 temp3 temp4];

%[h1,gap]=h1_without_gaps(h0,dh);
lh1=length(h1);
%h1=[h1,h1(lh1)+dh:dh:h1(lh1)+5*dh];

%dt=0.004;
pnorm=1;
option='Yilmaz ';
vel_min=1000;vel_max=3000;
np=120;dp=1/vel_min/np;
niter=10;
hyperparameter=1e-2;  %Use 1e-2 for pnorm=1; 1e+2 for pnorm=2
optionmute='n';
freqint=[2 256];
option_residuals='y';

[u,haxis,ttaxis,vr,ppaxis]=interp_taup1(x,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter,vel_max,optionmute,freqint,option_residuals);
figure,wigb(vr,1,ppaxis,ttaxis);
%[vel,pr]=alfa2vel(ppaxis,[0 dp]);
%pr=padzeros(pr,length(ppaxis));
%vel(1:10)=0;
%subplot(211),plot(vel);title('Vel');xlabel('#trace')
%subplot(212),plot(ppaxis);title('alfa');title('alfa');xlabel('#trace')




 
