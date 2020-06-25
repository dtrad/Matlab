% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=64;
coeff(1)=0.2;
coeff(2)=0.3;
coeff(3)=0.5;
dz(1)=100;
dz(2)=250;
dz(3)=300;
vel(1)=1500;
vel(2)=2000;
vel(3)=2300;
r0=-1;
dt=0.004;
t=0:NF-1;t=t*dt;
dh=10;
hh=0:NH-1;hh=hh*dh;
f=70;
nlayers=length(vel);
xa=dh*NH/2;
hh=hh-xa;
%xa=0;
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
%wav=[zeros(NF,NH/2) wav(:) zeros(NF,NH/2-1)];
S=fft(wav);
P1=zeros(NF,NH);
ss=zeros(1,NF);
ss(1)=1;
for nl=nlayers:-1:1
   r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft(r1);
   %R1=duplic(windowing(R1(1:NF/2,:)));
   %R1=duplic(R1(1:NF/2,:));

   if nl>1
   r2=zeros(NF,NH);coeff2=-coeff(nl-1);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
	elseif nl==1
   r2=zeros(NF,NH);coeff2=r0;r2(1,NH/2)=r0;R2=fft(r2);
	end

   for ii=1:NH
        	[W1]=propagator(dt,dh,NF,NH,vel(nl)/2,dz(nl),xa);
	end
   w1=ifft(W1);
   %w1=mix_fft(w1,wav);
   W1=fft(w1);
   P11=zeros(NF/2,NH*2-1);
  	for ii=4:NF/2-50;
   	P11(ii,:)=conv(W1(ii,:),(P1(ii,:)+R1(ii,:)));
   end
      
   %p11=ifft(duplic(windowing(P11(1:NF/2,:))));
   p11=ifft(duplic(P11(1:NF/2,:)));
   
   NHH=min(size(p11));
   p1=p11(1:NF,(NHH+1)/2-NH/2+1:(NHH+1)/2+NH/2);
   %p1=mix_fft(p1,wav);
   %p1=multiples2(p1,coeff2);
   P1=fft(p1);
   PF(:,:,nl)=P1;
end   
clear p   

for jj=1:nlayers   
	for ii=1:NH;
   	PF(:,ii,jj)=PF(:,ii,jj).*S(:);
	end
	p(:,:,jj)=ifft(PF(:,:,jj));
end

figure,
subplot(221),wigb(real(w1),1,hh,t);title('propagator');xlabel('offset');ylabel('time')
subplot(222),wigb(real(p(:,:,3)),1,hh,t);title('1 layer and halfspace');xlabel('offset');ylabel('time')
subplot(223),wigb(real(p(:,:,2)),1,hh,t);title('2 layer and halfspace');xlabel('offset');ylabel('time')
subplot(224),wigb(real(p(:,:,1)),1,hh,t);title('3 layer and halfspace');xlabel('offset');ylabel('time')


figure,
wigb(real(p(1:NF/2,:,1)),1,hh,t(1:NF/2));
title('Primaries- 3 layer model, (obtained with spatial convolution)')
xlabel('offset')
ylabel('time')




