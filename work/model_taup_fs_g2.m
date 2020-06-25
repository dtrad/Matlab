% Direct model of a seismogram in tau-p domain
% Using recursive formulas for Ri-1=f(Ri) the Ri traces are computed from the 
% halfspace (i=4) to surface (i=0) in the (w,p) domain. 
% It is assumed that all the energy comes from the surface so that R[nl](w,p)=0.
% Applying the recursive formula R[nl-1](w,p)=c3  
% R[i-1]=f(R[i])
% The Ro(tau,p) is obtained with the Inverse Fourier Transform
%
% Daniel Trad- Geop 520B- 1997.
% Based in INVERSION OF COMMON MIDPOINT SEISMIC DATA. (UFBa, Brazil, 1986)
clear;
close all;
% Known Parameters
optionxt='n';
option_noise='n';
optionzero='y';
optionsave='n';
option_direct_wave='n'

% velocities

V0=0;		
V(1)=1000;
V(2)=1500;
V(3)=2000;
%V(4)=3000;
%V(5)=4000;
% densities

r0=0;
r(1)=1000;
r(2)=1300;
r(3)=2000;
%r(4)=3000;
%r(5)=4000;

% two way travel distances;
% Three layer model with z1=300 m, z2=300 m and z3=400 m.

h(1)=400;
h(2)=600;
%h(3)=450;
%h(4)=600;

NF=512; % Number of frequencies
NP=80;  % Number of p traces.
NH=40;
dh=10;
dp=1/(NP*V(1)); % Delta p
dt=0.004;

nl=length(V);
w=frequency(dt,NF);
j=1:NP;
p=(j-1)*dp;

% First response R(nl)=0
ll=nl-1;
V1=V(ll);   
r1=r(ll);
V2=V(ll+1);   
r2=r(ll+1);
Y1=(sqrt(1-(p*V1).^2)/(r1*V1));
Y2=(sqrt(1-(p*V2).^2)/(r2*V2));
tt=ones(NF,1);
ttaxis=(0:NF-1)*dt;
ppaxis=p;
wav=rickerm(75,0.004);

if optionzero=='y'   % Generate zero offset output from Claerbout
   [x,cc]=trace(r,V,h,dt);x(NF)=0;x(1)=0;
   xc=convlim(x,wav,NF);
   X=[xc xc];
   figure,wigb(X,1,[1 2],ttaxis);
end   

c1=(Y1-Y2)./(Y1+Y2+eps);
mycomplexplot(c1,ppaxis);

c1=tt*c1;
R1=c1;

wp2taup2(R1,wav,ppaxis,ttaxis,'complex');
title('R initial = c1');

wp2taup2(R1,wav,ppaxis,ttaxis,'full');
title('R initial = c1 ');

for ll=nl-2:-1:1;
   R2=R1;
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
    
   dt2=(h2/V2)*((1-(p*V2).^2).^0.5);
   Y2=(sqrt(1-(p*V2).^2)/(r2*V2));
   if ll==0 
      c1=zeros(size(Y2));
      %c1=-0.8.*ones(size(Y2));
   elseif ll~=0
      c1=(Y1-Y2)./(Y1+Y2+eps);
   end;
   
   mycomplexplot(c1,ppaxis);
   c1=tt*c1;
   
   figure,
   subplot(211),plot(real(dt2));title('Imag dt2');
   subplot(212),plot(imag(dt2));title('Imag dt2');
   
   expterm=exp(-i*w*conj(dt2));
   
   R1=(c1+R2.*expterm)./(eps+1+c1.*R2.*expterm);clear expterm;
   
   %wp2taup2(R1,wav,ppaxis,ttaxis,'complex');
   
   %title(mytext);
  
   R1t=wp2taup2(R1,wav,ppaxis,ttaxis,'full');
   mytext=sprintf('Direct Model R%d(tau,p) %d layers',ll,nl-1);
   title(mytext);
end;
%[Q,D]=attenuation(ttaxis,ppaxis,V(1));
%R1t=R1t.*Q; % Attenuation


figure,
wigb(R1t,1,ppaxis,ttaxis);title('R1 Inside upper halfspace z=h(1)');

wavl=padzeros(wav,NF);
Aw=fft(wavl);

q=(1/V(1)^2-(p.^2)).^0.5;

w=w(:);
p=p(:).';
q=q(:).';
z=0;
z0=h(1)/2;
hs=100;
if option_direct_wave=='y'
P1=(Aw(:)*ones(1,NP)).*(exp(-i*w*q*z)+R1.*exp(-i*w*q*(2*z0-z)));
else 
P1=(Aw(:)*ones(1,NP)).*(R1.*exp(-i*w*q*(2*z0-z)));
end   

P1t=ifft(duplic(P1(1:NF/2,:)));

if option_noise=='y'
   PN=random('Normal',0,0.01,NF,NP);
   P1t=P1t+PN;
end   

figure,wigb(P1t,1,ppaxis,ttaxis);ylabel('time (sec)'),xlabel('p (s/m)');
mytext=sprintf('Direct Model, Pressure field P%d(t,x) %d layers',ll,nl-1);
title(mytext);


if optionxt=='y'
p=(0:NP-1)*dp;
h=(0:NH-1)*dh;
P1xt=taup2xt_inv(P1t(:,:),10,0.1,'linear',h,p(:),dt);figure,wigb(real(P1xt),1,h,ttaxis);
ylabel('time (sec)'),xlabel('offset (m)');
mytext=sprintf('Direct Model, Pressure field P%d(t,x) %d layers',ll,nl)
title(mytext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if optionsave=='y'
save c:\daniel\thesis\P1taup.mat P1t p h ttaxis
save c:\daniel\thesis\P1xt.mat P1xt p h ttaxis
end

z=10;
R1t=ifft(R1(1:NF/2,:));
R1t(NF/2-30,:)=0;

NDEC=2; 
NF=NDEC*NF;

temp=zeropadm(R1t,NF);
R1=fft(temp);clear temp;
w=frequency(dt,NF);
wavl=padzeros(wav,NF);
Aw=fft(wavl);
ttaxis=(0:NF-1)*dt;


GREEN_S=R1.*exp(-i*w*q*(2*z0-z-hs))./(1+R1.*exp(-2*i*z0*w*q)).*(exp(-2*i*hs*(w*q))-1).*(exp(-2*i*z*w*q)-1);
if option_direct_wave=='y'
   GREEN_INC=exp(-i*w*q*abs(z-hs))-exp(-i*w*q*(z+hs));
   GREEN=i./(2*w*q+eps).*(GREEN_INC+GREEN_S);
else
   GREEN=i./(2*w*q+eps).*(GREEN_S);
end
P0=(Aw(:)*ones(1,NP)).*GREEN;
P0t=ifft(duplic(P0(1:NF/2,:)));
figure,wigb(P0t(:,1:40),1,p(1:40),ttaxis);

if optionxt=='y'
p=(0:NP-1)*dp;
h=(0:NH-1)*dh;
P0xt=taup2xt_inv(P0t(:,1:40),10,0.1,'linear',h,p(1:40),dt);figure,wigb(real(P0xt),1,h,ttaxis);
ylabel('time (sec)'),xlabel('offset (m)');
mytext=sprintf('Direct Model, Pressure field P%d(t,x) %d layers',ll,nl)
title(mytext);
end
