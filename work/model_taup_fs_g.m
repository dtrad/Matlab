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
optionxt='n';
option_noise='n';
optionzero='y';
optionsave='n';
% velocities

V0=0;		
V(1)=1000;
V(2)=1300;
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

h(1)=200;
h(2)=300;
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

for ll=nl-2:-1:0;
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
   
   wp2taup2(R1,wav,ppaxis,ttaxis,'complex');
   mytext=sprintf('Direct Model R%d(tau,p) %d layers',ll,nl)
   title(mytext);
  
   R1t=wp2taup2(R1,wav,ppaxis,ttaxis,'full');
   title(mytext);
end;
%[Q,D]=attenuation(ttaxis,ppaxis,V(1));
%R1t=R1t.*Q; % Attenuation
Ro_nfs=R1t;
figure,wigb(Ro_nfs,1,ppaxis,ttaxis);title('Ro_nfs');

NDEC=1; % Number of surface multiples  to add
R1t=ifft(R1);
tempt=zeropadm(R1t,NDEC*NF);
temp=fft(tempt);clear tempt;
%Ro=temp.*(1+temp+temp.^2+temp.^3+temp.^4);

Ro=temp;
for ii=2:NDEC
   Ro=Ro+temp.^ii;
end


%Ro=R1./(1-R1);

%ttaxis=(0:4*NF-1).*dt;

Rot=wp2taup2(Ro,wav,ppaxis,ttaxis,'complex');
title(mytext);

Rot=wp2taup2(Ro,wav,ppaxis,ttaxis,'full');
title(mytext);

%Rot=normalize(Rot);
if option_noise=='y'
   RN=random('Normal',0,0.01,NF,NP);
   Rot=Rot+RN;
end   


if optionxt=='y'
p=(0:NP-1)*dp;
h=(0:NH-1)*dh;
Ro=taup2xt_inv(Rot(:,:),10,0.1,'linear',h,p(:),dt);figure,wigb(real(Ro),1,h,ttaxis);
ylabel('time (sec)'),xlabel('offset (m)');
mytext=sprintf('Direct Model R%d(t,x) %d layers',ll,nl)
title(mytext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if optionsave=='y'
save c:\daniel\thesis\Rotaup.mat Rot Ro_nfs p h ttaxis
save c:\daniel\thesis\Roxt.mat Ro p h ttaxis
end

wavl=padzeros(wav,NF);
Aw=fft(wavl);

Vo=500;
q=(1/Vo^2-(p.^2)).^0.5;

w=w(:);
p=p(:).';
q=q(:).';
z=0;
z0=100;
hs=10;
Ro=temp;clear temp;
%P0=(Aw(:)*ones(1,NP)).*(0*exp(-i*w*q*z)+Ro.*exp(-i*w*q*(2*z0-z)));

GREEN_INC=exp(-i*w*q*abs(z-hs))-exp(-i*w*q*(z+hs));
GREEN_S=Ro.*exp(-i*w*q*(2*z0-z-hs))./(1+Ro.*exp(-2*i*z0*w*q)).*(1-exp(-2*i*hs*(w*q))).*(1-exp(-2*i*z*w*q));

GREEN=i/(2*w*q).*(GREEN_INC+GREEN_S);
P0=(Aw(:)*ones(1,NP)).*GREEN;


P0t=ifft(duplic(P0(1:NF/2,:)));

P0xt=taup2xt_inv(P0t(:,:),10,0.1,'linear',h,p(:),dt);figure,wigb(real(P0xt),1,h,ttaxis);
