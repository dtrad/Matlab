% Seismology: Lab #1
% Daniel Trad
clear
alpha=6000;
beta=3500;
rho=2700;
mu_s=0.2;

mu=beta^2*rho;

lambda=alpha^2*rho-2*mu;

el(1,1)=-0.26e-6;
el(2,2)=0.92e-6;
el(1,2)=-0.69e-6;
el(2,1)=el(1,2);

taul=stress(lambda,mu,el);

[V,D]=eig(el);

ey(1,1)=0.101e-6;
ey(2,2)=-0.020e-6;
ey(1,2)=0.05e-6;
ey(2,1)=ey(1,2);

tauy=stress(lambda,mu,ey);
tau1000y=1000*tauy;

tau=tau1000y;

if (0) 
tau(1,1)=9.196e6;
tau(1,2)=3.31e5;
tau(2,2)=1.192e6;
tau(2,1)=3.31e5;
end
ang=[0:10:170];
nz=length(ang);

for iz=1:nz

  theta=ang(iz)*pi/180;
  n=[cos(theta);-sin(theta)];
  s=[sin(theta);cos(theta)];
  traction=tau*n;
  
  
  tn(iz)=traction'*n;
  ts(iz)=traction'*s;
  dcff(iz)=abs(ts(iz))+mu_s*(tn(iz));
  
end,  
figure;
subplot(311); plot(ang,tn,'o');title('Seislab');
subplot(312); plot(ang,ts,'o');
subplot(313); plot(ang,dcff,'o');


