% Load data generated with synth0. xxc contains dh,dt and h_near
%load c:\daniel\thesis\xxc;
z=seis_shape(xxc(1:256,1:40));clear data;
z=normalize(z);
%dh=25;
select=ones(1,40);
x=shrinkt2(z,select);clear z;
temp=0:length(select)-1;temp=temp*dh;
h0=shrinkt2(temp,select);clear temp;
[nt,nh]=size(x);
h1=h0;
dt=0.004;
pnorm=1;
option='Yilmaz ';
vel_min=2000;vel_max=5000;
np=50;niter=5;
hyperparameter=1e-2;  %Use 1e-2 for pnorm=1; 1e+2 for pnorm=2
optionmute='n';
freqint=[2 127];
option_residuals='y';
[u,haxis,ttaxis,vr,ppaxis]=interp_taup1(x,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter,vel_max,optionmute,freqint,option_residuals);
figure,wigb(vr,1,ppaxis,ttaxis);



 