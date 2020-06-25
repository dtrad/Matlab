% Load data generated with synth0. xxc contains dh,dt and h_near
load c:\daniel\thesis\xxc;

z=seis_shape(xxc(1:256,1:50));clear data;
z=normalize(z);
%dh=25;
select=[ones(1,10) ones(1,10) zeros(1,0) ones(1,30) zeros(1,0)];
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

[h1,gap]=h1_without_gaps(h0,dh);
lh1=length(h1);
h1=[h1,h1(lh1)+dh:dh:h1(lh1)+0*dh];

%dt=0.004;
pnorm=1;
option='Hampson';
vel_min=2500;vel_max=3300;
np=80;dp=1/vel_min/np;
niter=10;
hyperparameter=1e-2;  %Use 1e-2 for pnorm=1; 1e+2 for pnorm=2
optionmute='y';
freqint=[2 127];
option_residuals='n';

[u,haxis,ttaxis,vr,ppaxis]=interp_taup1(x,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter,vel_max,optionmute,freqint,option_residuals);
figure,wigb(vr,1,ppaxis,ttaxis);
%[vel,pr]=alfa2vel(ppaxis,[0 dp]);
%pr=padzeros(pr,length(ppaxis));
%vel(1:10)=0;
%subplot(211),plot(vel);title('Vel');xlabel('#trace')
%subplot(212),plot(ppaxis);title('alfa');title('alfa');xlabel('#trace')
plotwig_mult(x,vrorig,u,vr,ttaxis(1:256),ppaxis,h0,h1);

 