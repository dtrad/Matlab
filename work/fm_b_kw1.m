% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=64;
coeff(1)=0.2;
coeff(2)=0.3;
coeff(3)=-0.5;
dz(1)=300;
dz(2)=550;
dz(3)=300;
vel(1)=1000;
vel(2)=2000;
vel(3)=2300;

dt=0.004;
dh=10;
hh=0:NH-1;hh=hh*dh;
f=70;
nlayers=length(vel);
xa=dh*NH/2;
%xa=0;
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
%wav=[zeros(NF,NH/2) wav(:) zeros(NF,NH/2-1)];
S=fft(wav);S=S(:)*ones(1,NH);
P1=zeros(NF,NH);

for nl=nlayers:-1:1
   r1=zeros(NF,NH);r1(1,1)=coeff(nl);R1=fft2(r1);
   for ii=1:NH
        	[W1]=propagator(dt,dh,NF,NH,vel(nl),dz(nl),xa);
	end
   w1=ifft(W1);
   W1=fft2(w1);
   P1=[P1 zeros(size(P1))];
   P1=[P1;zeros(size(P1))];
   R1=[R1 zeros(size(R1))];
   R1=[R1;zeros(size(R1))];
   W1=[W1 zeros(size(W1))];
   W1=[W1;zeros(size(W1))];
   

   P1=W1.*(P1+R1);
   p1=ifft2(P1);
   
   p1=p1(1:NF,2*NH/4+1:2*3*NH/4);
   P1=fft2(p1);
   figure,wigb(p1);
   PF(:,:,nl)=P1;
   
end   
   
for jj=1:nlayers   
  	PF(:,:,jj)=PF(:,:,jj).*S;
	p(:,:,jj)=ifft2(PF(:,:,jj));
end

figure,
subplot(221),wigb(real(w1));
subplot(222),wigb(real(p(:,:,1)));
subplot(223),wigb(real(p(:,:,2)));
subplot(224),wigb(real(p(:,:,3)));

figure,
wigb(real(p(1:NF/2,:,1)));



