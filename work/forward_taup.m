% Direct model of a seismogram in tau-p domain
% Using recursive formulas for Ri-1=f(Ri) the Ri traces are computed from the 
% halfspace (i=4) to surface (i=0) in the (w,p) domain. 
% It is assumed that all the energy comes from the surface so that R4(w,p)=0.
% Applying the recursvive formula R3(w,p)=c3  
% R2=f(R3), R1=f(R2) and Ro=f(R1)
% The Ro(tau,p) is obtained with the Inverse Fourier Transform
%
% Daniel Trad- Geop 520B- 1997.
% Based in Dimas Ferreyra's thesis. (UFBa, Brazil, 1986)
clear;
close all;
% Known Parameters

% velocities

V0=0;		
V1=1000;
V2=2200;
V3=4500;
V4=6300;

% densities

r0=0;
r1=1100;
r2=2200;
r3=3300;
r4=4400;

% two way travel distances;
% Three layer model with z1=300 m, z2=300 m and z3=400 m.

h1=200;
h2=300;
h3=100;

% Nyquist Freq. = 512 Hertz
% Sample freq.=1024 Hertz
% dt = 1/1024 = 0.001 ;
% densities in arbitrary units.
% delta p = 1/20000 = 5.e-5; 

NF=512; % Number of frequencies
NP=60;  % Number of p traces.
NH=40;
dh=25;
dp=1/(NP*V1); % Delta p
dt=1/1024;FN=512;
% because we need one p at each time the functions are not saved

for j=1:NP;
  p=(j-1)*dp;

% Admittances in function of p.
% Yo=0;  % (free surface) 
% Yo=Y1; % (no free surface)

  Y1=sqrt(1-(p*V1)^2)/(r1*V1);
  Y2=sqrt(1-(p*V2)^2)/(r2*V2);
  Y3=sqrt(1-(p*V3)^2)/(r3*V3);
  Y4=sqrt(1-(p*V4)^2)/(r4*V4);

% Reflection coefficients.
% co=(Yo-Y1)/(Yo+Y1); model with free surface.
% co=0; model without free surface. (Yo=Y1).

 co=0;  
 c1=(Y1-Y2)/(Y1+Y2);
 c2=(Y2-Y3)/(Y2+Y3);
 c3=(Y3-Y4)/(Y3+Y4);

 dt1=(h1/V1)*((1-(p*V1)^2)^0.5);
 dt2=(h2/V2)*((1-(p*V2)^2)^0.5);
 dt3=(h3/V3)*((1-(p*V3)^2)^0.5);

 m=1:NF;
 w=2*pi*(m'-1); % the Nyquist Freq. is 512 hertz--> dt=1/1024
 R4=0;
 R3=c3;
 R2=(c2+R3.*exp(i*w*dt3))./(1+c2.*R3.*exp(i*w*dt3));
 R1=(c1+R2.*exp(i*w*dt2))./(1+c1.*R2.*exp(i*w*dt2));
 Ro(m,j)=(co+R1.*exp(i*w*dt1))./(1+co.*R1.*exp(i*w*dt1));
 

% The Ro(w,p) values are arranged to form the FFT of a real trace.

m=1:NF;
end;

Rop=duplic(Ro);
Rot=ifft(Rop);
Rot=reverse(Rot);

figure,
wigb(Rot);
xlabel('time sec'),ylabel('# p-trace');
title('direct model Ro(tau,p)- three layers');

p=(0:NP-1)*dp;
h=(0:NH-1)*dh;

Ro=taup2xt_inv(Rot,[],[],'linear',h,p,dt);figure,wigb(real(Ro));

%Ro=taup2xt(Rot,'linear',h,p,dt);
%figure,wigb(Ro);

