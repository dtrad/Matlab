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
vel(2)=1000;
dt=0.004;
dh=10;
f=70;
nlayers=length(vel);
t=0:NF-1;t=t*dt;
h=0:NH-1;h=h*dh;

wav=rickerm(f,dt);
wav=padzeros(wav,NF);
wav=[wav(:) zeros(NF,NH-1)];
S=fft2(wav);
figure,wigb(wav);title('source');
P=zeros(NF,NH);

for nl=nlayers:-1:1
   
   % Construct the propagator
   [kz,kx,w,aa]=obliquity(dt,dh,NF,NH,vel(nl));
	kz=kz(1:NF,1:NH);
	W1=exp(-i*kz*2*dz(nl));
   
   % Delete the direct wave
   %kk=(abs(kz)>1e-10);
   %W1=W1.*kk;
   W1=W1.*aa;
   W1=zeros_2d(W1); % Clean the right dipping events
   w1=ifft2(W1);
   %w1(1:25,:)=0;
   w1(:,NH/2+1:NH)=0;
   %w1(NF/2+2:NF,:)=0;
   figure(1+nl),
   subplot(223),
   ww=fftshift(w);kxx=fftshift(kx);ww(1)=ww(2);kxx(1)=kxx(2);
   wigb(fftshift(W1),1,kxx,ww);title('propagator in kx-w');
   xlabel('kx')
   ylabel('freq')

   subplot(224),
   wigb(w1,1,h,t);title('propagator in x-t');
   xlabel('offset')
   ylabel('time')

   W1=fft2(w1);
%  P=(W1.*(P+R1).*W1);

	r1=zeros(NF,NH);r1(1,1)=coeff(nl);R1=fft2(r1);

	P=(W1.*(P+R1));
   PM=P./(1+R1.*P);
   if nl==1 P=P.*S;end
   p=ifft2(P);
   pm=ifft2(PM);
   %p=ifft2(duplic2d(P));
   figure(1+nl),
   subplot(221),
   wigb(p,1,h,t);title('synthetic data: primaries') 
   xlabel('offset')
   ylabel('time')

   subplot(222)
   wigb(pm,1,h,t);
   title('synthetic data: primaries+multiples') 
   xlabel('offset')
   ylabel('time')

end


