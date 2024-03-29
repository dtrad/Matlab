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
h(2)=400;
%h(3)=450;
%h(4)=600;

NF=512; % Number of frequencies
NP=80;  % Number of p traces.
NH=40;
dh=25;
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
pmax=1/V2;
ppmax=t2index(pmax,dp);
ppmax

Y1=(sqrt(1-(p(1:ppmax)*V1).^2)/(r1*V1));
Y2=(sqrt(1-(p(1:ppmax)*V2).^2)/(r2*V2));
tt=ones(NF,1);
ttaxis=(0:NF-1)*dt;
ppaxis=p;

c1=(Y1-Y2)./(Y1+Y2+eps);
c1(ppmax+1:NP)=zeros(1,NP-ppmax);
mycomplexplot(c1,ppaxis);

c1=tt*c1;


R1=c1;
wav=rickerm(50,0.004);
R1t=ifft(c1);

for ii=1:NP,R1t(:,ii)=convlim(R1t(:,ii),wav,NF);end
figure,wigb(R1t,1,ppaxis,ttaxis);title('R initial = c1');
xlabel('p (s/m)');ylabel('t (sec)');

wp2taup2(R1,wav,ppaxis,ttaxis,'full');
title('R initial = c1');

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
 
   if ll~=0 pmax=1/V2;ppmax=t2index(pmax,dp);elseif ll==0 ppmax=NP;end
   ppmax
   dt2=(h2/V2)*((1-(p*V2).^2).^0.5);
   Y2=(sqrt(1-(p*V2).^2)/(r2*V2));
   if ll==0 
      c1=zeros(size(Y2));
      %c1=-0.8.*ones(size(Y2));
   elseif ll~=0
      c1=(Y1-Y2)./(Y1+Y2+eps);
   end;
   c1=c1;
   mycomplexplot(c1,ppaxis);
   c1(ppmax+1:NP)=zeros(1,NP-ppmax);
   
   c1=tt*c1;
   %c1=vector2matrix(c1,dt2,dt,NF);
   %figure,wigb(real(c1));
   %c1=fft(c1);
   figure,
   subplot(211),plot(real(dt2));title('Imag dt2');
   subplot(212),plot(imag(dt2));title('Imag dt2');
   
   temparg=w*conj(dt2);
   temparg(NF/2+2:NF,:)=conj(temparg(NF/2+2:NF,:));
   tempexp=exp(-i*temparg);
   R1=(c1+R2.*tempexp)./(eps+1+c1.*R2.*tempexp);
   
   R1t=ifft(R1);
   for ii=1:NP,R1t(:,ii)=convlim(R1t(:,ii),wav,NF);end
  %RN=random('Normal',0,0.01,NF,NP);
  %R1t=R1t+RN;
   NNP=NP;
   figure,wigb(R1t(:,1:NNP),1,p(1:NNP),ttaxis);
   ylabel('time (sec)'),xlabel('p (s/m)');
   mytext=sprintf('Direct Model R%d(tau,p) %d layers',ll,nl)
   title(mytext);
   wp2taup2(R1,wav,ppaxis,ttaxis,'full');
   title(mytext);

end;
temp=R1;
%Ro=R1.*(1+R1.*exp(-i*w*conj(dt2))-R1.^2.*exp(-i*w*2*conj(dt2)));
Ro=R1.*(1+temp-temp.^2+temp.^3);
Rot=ifft(Ro);
for ii=1:NP,Rot(:,ii)=convlim(Rot(:,ii),wav,NF);end

figure,wigb(Rot);
%Rot=R1t;
Rot=wp2taup2(Ro,wav,ppaxis,ttaxis,'full');
mytext=sprintf('Direct Model R%d(tau,p) %d layers',ll,nl)
title(mytext);
save c:\daniel\synth\Rot.mat Rot

if optionxt=='y'
p=(0:NP-1)*dp;
h=(0:NH-1)*dh;
Ro=taup2xt_inv(Rot,[],[],'linear',h,p,dt);figure,wigb(real(Ro),1,h,ttaxis);
ylabel('time (sec)'),xlabel('offset (m)');
mytext=sprintf('Direct Model R0(t,x) %d layers',nl)
title(mytext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
