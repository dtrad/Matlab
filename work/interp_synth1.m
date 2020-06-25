% Load data generated with synth0. xxc contains dh,dt and h_near
load /home/dtrad/matlab/xxc;
x=seis_shape(xxc(1:256,1:64));
h=h(1:64);

clear xxc;
x=normalize(x);

[nt,nh]=size(x);

dh=25;

h0=h;
h1=h(1):dh:h(nh);
dt=0.004;
pnorm=1;
option='Hampson';
vel_min=1300;vel_max=2500;
np=40;dp=1/vel_min/np;
niter=5;
hyperparameter=1e-2;  %Use 1e-2 for pnorm=1; 1e+2 for pnorm=2
optionmute='n';
freqint=[2 127];
option_residuals='n';

[u,haxis,ttaxis,vr,ppaxis]=radon1(x,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter,vel_max,optionmute,freqint,option_residuals);
figure,wigb(vr,1,ppaxis,ttaxis);
%[vel,pr]=alfa2vel(ppaxis,[0 dp]);
%pr=padzeros(pr,length(ppaxis));
%vel(1:10)=0;
%subplot(211),plot(vel);title('Vel');xlabel('#trace')
%subplot(212),plot(ppaxis);title('alfa');title('alfa');xlabel('#trace')




 





