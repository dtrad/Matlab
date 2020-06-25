% Direct model of a seismogram in tau-p domain
% Using recursive formulas for Ri-1=f(Ri) the Ri traces are computed from the 
% halfspace (i=4) to surface (i=0) in the (w,p) domain. 
% It is assumed that all the energy comes from the surface so that R[nl](w,p)=0.
% Applying the recursvive formula R[nl-1](w,p)=c3  
% R[i-1]=f(R[i])
% The Ro(tau,p) is obtained with the Inverse Fourier Transform
%
% Daniel Trad- Geop 520B- 1997.
% Based in Dimas Ferreyra's thesis. (UFBa, Brazil, 1986)
clear;
close all;
% Known Parameters

% velocities

V0=0;		
V(1)=1000;
V(2)=1500;
V(3)=2000;
%V(4)=6300;
% densities

r0=0;
r(1)=1000;
r(2)=1500;
r(3)=2000;
%r(4)=4400;

% two way travel distances;
% Three layer model with z1=300 m, z2=300 m and z3=400 m.

h(1)=200;
h(2)=300;
%h(3)=100;

NF=512; % Number of frequencies
NP=80;  % Number of p traces.
NH=40;
dh=25;
dp=1/(NP*V(1)); % Delta p
dt=0.004;

nl=length(V);
m=1:NF;
w=2*pi*(m'-1);
j=1:NP;
p=(j-1)*dp;

% First response R(nl)=0
ll=nl-1;
V1=V(ll);   
r1=r(ll);
V2=V(ll+1);   
r2=r(ll+1);
Y1=sqrt(1-(p*V1).^2)/(r1*V1);
Y2=sqrt(1-(p*V2).^2)/(r2*V2);
tt=ones(NF,1);
c1=(Y1-Y2)./(Y1+Y2);c1=tt*c1;
R1=c1;

for ll=nl-2:-1:0;
 R2=R1;
 if ll~=0 V1=V(ll);elseif ll==0 V1=0;end
 if ll~=0 r1=r(ll);elseif ll==0 r1=0;end
   
 V2=V(ll+1);   
 r2=r(ll+1);
 h2=h(ll+1);
  
 if ll~=0 Y1=sqrt(1-(p*V1).^2)/(r1*V1);elseif ll==0 Y1=zeros(size(p));end
 Y2=sqrt(1-(p*V2).^2)/(r2*V2);

 tt=ones(NF,1);
 co=0;  
 c1=(Y1-Y2)./(Y1+Y2);if ll==0 c1=zeros(size(Y2));end;c1=tt*c1;
 dt2=(h2/V2)*((1-(p*V2).^2).^0.5);
 R1=(c1+R2.*exp(i*w*dt2))./(1+c1.*R2.*exp(i*w*dt2));

end; 

Ro=R1;
Rop=duplic(Ro);
Rot=ifft(Rop);
Rot=reverse(Rot);
wav=rickerm(30,0.004);
for ii=1:NP,Rot(:,ii)=convlim(Rot(:,ii),wav,2*NF);end
RN=random('Normal',0,0.01,2*NF,NP);
Rot=Rot+RN;

figure,
wigb(Rot);
xlabel('time sec'),ylabel('# p-trace');
title('direct model Ro(tau,p)- three layers');

p=(0:NP-1)*dp;
h=(0:NH-1)*dh;

Ro=taup2xt_inv(Rot,[],[],'linear',h,p,dt);figure,wigb(real(Ro));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
