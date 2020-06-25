load c:\daniel\thesis\P0xt;
z=seis_shape(P0xt(1:512,:));clear P0xt;
z=normalize(z);

[nt,nh]=size(z);

%dh=25;
%temp1=(0:9)*dh;
%temp2=(10:11)*dh; % Gap
%temp3=(12:29)*dh;

%h0=[temp1 temp2 temp3];
%h1=[temp1 temp2 temp3];

h0=h;
h1=h;


%select=[ones(1,10) ones(1,2) ones(1,18)];
%x=shrinkt2(z,select);
x=z;clear z;

dt=0.004;
pnorm=1;
option='Hampson';
vel_min=500;
np=60;
niter=5;
hyperparameter=1e-5;

[u,haxis,ttaxis,vr,ppaxis]=interp_taup(x,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter);
figure,wigb(vr,1,ppaxis,ttaxis);
