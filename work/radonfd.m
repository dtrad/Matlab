close all
clear
synthrad;
vmin=1000;
vmax=3500;
perc=50;
method=3;
reorth=1;

iterend=5;
eps1=1e-3;
eps2=1e-3;
inneriter=2;
step=.9;

[nt np]=size(v);
[nt nh]=size(d);
noise=0.5*perc/100*(rand(size(d))-0.5);

d=d+noise;
figure,wigb(d,1,h,t);

%dr=radonop(v,t,h,p);
[vr,dr]=radon0(d,h,np,vmin,vmax,dt,method,[],eps1,eps2,.9);
figure,wigb(vr,1,p,t);  
  
for i=1:iterend;
[vr,dr]=radon0(d,h,np,vmin,vmax,dt,method,vr,eps1,eps2,.9);
figure,wigb(vr,1,p,t);  
end


