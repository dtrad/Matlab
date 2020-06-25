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
vel(1)=500;
vel(2)=1000;
dt=0.004;
dh=25;
f=70;
nlayers=length(vel);
xa=dh*NH/2;
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
wav=[wav(:) zeros(NF,NH-1)];
S=fft(wav);
P=zeros(NF,NH);

for nl=nlayers:-1:2
   %[kz,kx,w]=obliquity(dt,dh,NF,NH,vel(nl));
	%kz=kz(1:NF,1:NH);
	%[kz,theta,k]=kz_fxz(dz,dh,NH,w,vel);
	r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft(r1);
  	W1=propagator(dt,dh,NF,NH,vel(nl),dz(nl),xa);

	%W1=exp(-i*kz*dz(nl));
   w1=ifft(W1);
   w1(NF/2+1:NF,:)=0;
   figure,wigb(w1);
   W1=ifft(w1);
   %w1=[w1 zeros(size(w1))];
   %w1=[w1;zeros(size(w1))];
   
   %P2=fft2(w1).*fft2(r1).*fft2(w1);
   %p2=ifft2(P2);
   %W1=fft2(w1);
   for ii=1:NF;
   %P=(W1(1:NF,1:NH)*((P+R1).'));%.'*W1);
   P1(ii,:)=conv(W1(ii,:),R1(ii,:));
   end
   %P=P./(1+0*R1.*P);
   %if nl==1 P=P.*S;end
	p1=ifft(P1);
   figure,
   wigb(p1);
end


