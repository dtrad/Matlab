%mex -g  c:\daniel\cpp\radon_cary2.cpp

%load xxc;
%u=xx(1:256,1:30);
%clear xxc;clear xx;
v=4000;
u=[zeros(50,30);ones(1,30);zeros(77,30)];
%u(:,10)=1;
dt=0.004;
dh=25;
np=25;
dp=1/v/(np-1);
[v]=radonop1(u,dt,dh,dp,np);
figure,wigb(v);
