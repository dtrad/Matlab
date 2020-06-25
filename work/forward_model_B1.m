% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=64;
coeff(1)=0.55;
coeff(2)=0.55;
dz(1)=100;
dz(2)=250;
vel(1)=1000;
vel(2)=2000;
dt=0.004;
dh=25;
f=70;
nlayers=length(vel);

wav=rickerm(f,dt);
wav=padzeros(wav,NF);
wav=[wav(:) zeros(NF,NH-1)];
S=fft2(wav);
P=zeros(NF,NH);

for nl=1:nlayers
	[kz,kx,w]=obliquity(dt,dh,NF,NH,vel(nl));
	kz=kz(1:NF,1:NH);
	%[kz,theta,k]=kz_fxz(dz,dh,NH,w,vel);
	r1=zeros(NF,NH);r1(1,1)=coeff(nl);R1=fft2(r1);

	W1=exp(-i*kz*dz(nl));
   w1=ifft2(W1);
   w1(:,NH/2+1:NH)=0;
   w1(NF/2+1:NF)=0;
   W1=fft2(w1);
   P=P+(W1.*R1.*W1).*S;
   P=P./(1+R1.*P);
   S=P;

	p=ifft2(duplic2d(P));
	figure,wigb(p);
end


