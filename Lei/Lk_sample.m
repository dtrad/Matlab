clear;clc
dz=8;
dx=8;
nz=5;
nx=10;

h=[dz dx];
n=[nz nx];
parameterization=2;
rho=1.+zeros(5,10);
density_inv(:,1)=1./rho(:);
density=1./density_inv;
f=10;
model=zeros(5,10);
model(1:2,:)=1./2000;
model(3:5,:)=1./4000;
%imagesc(model);
sqmodel=model.^2;
m(:,1)=sqmodel(:);

Lk = Helmoholtz2D_vd_parameterization2(f,m,density_inv,h,n,parameterization,model);

 
