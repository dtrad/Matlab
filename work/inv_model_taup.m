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
V(2)=1300;
V(3)=2000;
V(4)=3000;
V(5)=4000;
% densities

r0=0;
r(1)=1000;
r(2)=1300;
r(3)=2000;
r(4)=3000;
r(5)=4000;
% two way travel distances;
% Three layer model with z1=300 m, z2=300 m and z3=400 m.

h(1)=200;
h(2)=300;
h(3)=450;
h(4)=600;
load c:\daniel\synth\Rot;
Rot=seis_shape(Rot);

[NF,NP]=size(Rot);
% NF Number of frequencies and time samples
% NP Number of p traces.
NH=40;
dh=25;
dp=1/(NP*V(1)); % Delta p
dt=0.004;

nl=length(V);
w=frequency(dt,NF);
j=1:NP;
p=(j-1)*dp;
pp=1:NP;
tt=ones(NF,1);
ttaxis=(0:NF-1)*dt;
ppaxis=p;

wav=rickerm(50,dt);
wavf=fft(wav);

lw=length(wav);

figure,wigb(Rot,1,ppaxis,ttaxis);ylabel('time sec'),xlabel('p (s/m)');
title('direct model Ro(tau,p)- five layers');

R2=fft(Rot);  % Data
for ll=0:nl-2;
   R1=R2;%R1=normalize(R1);
   R1=detrend(R1);
   if ll==0 
      V1=0;
      r1=0;
      Y1=zeros(size(p));
   elseif ll~=0
      V1=V(ll);
      r1=r(ll);
      Y1=(sqrt(1-(p*V1).^2)/(r1*V1));
   end
    
   V2=V(ll+1);   
 	r2=r(ll+1);
 	h2=h(ll+1);
  	Y2=(sqrt(1-(p*V2).^2)/(r2*V2));
   pmax=1/V(ll+1);  % Critic Angle
   ppmax=t2index(pmax,dp)-1;  % Index before critic angle
     
   if ll==0
      c1=zeros(size(Y2(1:ppmax)));
   elseif ll~=0
      c1theor=(Y1-Y2)./(Y1+Y2+eps);  % Reflection coeff. computed from the real parameters
      %c1=c2vector;                  % Reflection coeff. computed from the data 
      c1=c1theor(1:ppmax);
     	figure
   	subplot(211),plot(real(c1)),title('c1 from data');
   	subplot(212),plot(real(c1theor)),title('c1 theoretical');
   end;
   
   clear c1temp;

   c1=tt*c1;
   dt2=((h2/V2)*((1-(p(1:ppmax)*V2).^2).^0.5));  % Delta tau between Zn+1-Zn
   
   tempM=c1;
   tempMt=ifft(tempM);
   for ii=1:length(tempMt(1,:));tempMt(:,ii)=convlim(tempMt(:,ii),wav,NF);end;
   
   %figure,
   %subplot(211);plot(real(tempMt(:,1)));title('c1*wav');    
   %subplot(212);plot(imag(tempMt(:,1)));title('Imag c1*wav'); 
   %figure,wigb(tempMt);title('c1');
   
   
   c1temp=fft(tempMt);
   temp=c1.*R1(:,1:ppmax)./2;
   
   %R2=(R1(:,1:ppmax)-c1temp).*(1+temp-temp.^2+temp.^3).*exp(i*w*(dt2));
   MULT=(R1(:,1:ppmax)-c1temp).*(temp-temp.^2+temp.^3-temp.^4).*exp(i*w*(dt2));
   R2=(R1(:,1:ppmax)-c1temp)./(1-c1.*R1(:,1:ppmax)./2).*exp(i*w*(dt2));

   mult=ifft(MULT);
   
   figure,
   subplot(211);plot(real(mult(:,1)));title('real(multiple) order 1 p=0');    
   subplot(212);plot(imag(mult(:,1)));title('Imag(multiple) order 1 p=0'); 

   figure,wigb(mult);title('Mutiple of order 1');
   
%   XX=(1-c1.^2).*R1./(1-c1.*R1).*exp(i*w*(dt2));
%   R2=XX-c1;
   c2=ifft(R2);
   
   figure,
   subplot(211),plot(real(c2(:,1)));title('c2=ifft(R2)');
   subplot(212);plot(dt2(1:ppmax));title('dt2')
   c2vector=c2(round(lw/2),pp(1:ppmax));
   R2t=ifft(R2); 
   NNP=ppmax;
   
   R2t=normalize(R2t);
   figure,wigb(R2t(:,1:NNP),1,p(1:NNP),ttaxis);
   ylabel('time (sec)'),xlabel('p (s/m)');
   mytext=sprintf('Inverse Model R%d(tau,p) %d layers',ll+1,nl)
   title(mytext);
end; 

%p=(0:NP-1)*dp;
%h=(0:NH-1)*dh;
%Ro=taup2xt_inv(R2t,[],[],'linear',h,p,dt);figure,wigb(real(Ro));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
