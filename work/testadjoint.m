% Script to calculate the adjoint of a matrix
%load model.mat 
close all
np=64;
nt=512;
dt=0.004;

model=zeros(nt,np);model(nt/2,np/2)=1;

w=ricker(dt,35);
for i=1:np;mm(:,i)=conv(w,model(:,i));end;

model=mm(1:nt,1:np);


[nt np]=size(model);
rtmethod='LRT'
%define  axes

vmin=1000;
vmax=4000;


hmin=-2500;
hmax=-hmin;
nh=256;

h=hmin:(hmax-hmin)/(nh-1):hmax;

t=0:nt-1;t=t*dt;

p=radonaxis2(np,vmin,vmax,h,rtmethod);

if (1)
  [FF3,F3]=smearingfilt(t,h,p,rtmethod);
  save smearingmatrix.mat FF3 F3
else 
  load smearingmatrix.mat
end

[madj]=applysmearing(model,t,p,FF3);


figure(1); wigb(model);title('model');
figure(2); wigb(madj);title('madj');

[dp1]=applyrtop(model,t,h,F3);
figure(3); wigb(dp1);title('dp1');



[dp2]=applyrtop(madj,t,h,F3);
figure(4); wigb(dp2);title('dp2');


