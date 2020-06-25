% Script to calculate the adjoint of a matrix
%load model.mat 
close all
clear
np=128;
nt=128;
nh=128;
dt=0.004;
eps1=0.1;
perc=0.005;


model=zeros(nt,np);model(nt/2,np/2)=1;

if (0)
  w=ricker(dt,35);
  for i=1:np;mm(:,i)=conv(w,model(:,i));end;
  model=mm(1:nt,1:np);
end


[nt np]=size(model);
rtmethod='LRT'
%define  axes

vmin=1000;
vmax=4000;


hmin=-500;
hmax=-hmin;


h=hmin:(hmax-hmin)/(nh-1):hmax;

t=0:nt-1;t=t*dt;

p=radonaxis2(np,vmin,vmax,h,rtmethod);

if (1)
  [FF3,F3]=smearingfilt(t,h,p,rtmethod);
  save smearingmatrix.mat FF3 F3
else 
  load smearingmatrix.mat
end

[x]=applyrtop(model,t,h,F3);

noise=randn(nt,nh)*perc;
y=x+noise;
figure(1); wigb(y);title('y');

[madj]=applyadj(y,t,p,F3);

[mest]=applydesmearing(madj,t,p,FF3,eps1);
mest=real(mest);

figure(2),wigb(madj);title('madj')
figure(3),wigb(mest);title('mest')

%%%%%%%%% do the wavelet denoising


%^T = 0.0016;                     % noise level

T = 0.001;                     % noise level

[wmest]=applywt2(mest,nt,np,T);
figure(4),clf;wigb(wmest);title('wmest')

% data prediction

[dmest]=applyrtop(mest,t,h,F3);
[dwmest]=applyrtop(wmest,t,h,F3);

figure(5); wigb(dmest);title('mest');
figure(6); wigb(dwmest);title('dwmest');



%[madj]=applysmearing(model,t,p,FF3);

%[mp]=apply(madj,t,h,FF3);

return;


figure(1); wigb(model);title('model');
figure(2); wigb(madj);title('madj');




[dp2]=applyrtop(madj,t,h,F3);
figure(4); wigb(dp2);title('dp2');







