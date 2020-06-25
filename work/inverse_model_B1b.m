% Inverse Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=64;
r0=-0.6; % Surface reflection coefficient
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

r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft2(r1);
if nl>1
	r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft2(r2);
elseif nl==1
   r2=zeros(NF,NH);r2(1,NH/2)=r0;R2=fft2(r2);
end;
	P=(W1.*(P+R1));
   PM=P./(1-R2.*P);
   if nl==1 P=P.*S;P0=PM;PM=PM.*S;end
   p=ifft2(P);
   pm=ifft2(PM);
   %p=ifft2(duplic2d(P));
   figure(1+nl),
   subplot(221),
   wigb(p,1,h,t);
   mytext=sprintf('Primary from layer %d and P+M from below',nl); 
   title(mytext) 
   xlabel('offset')
   ylabel('time')

   subplot(222)
   wigb(pm,1,h,t);
   mytext=sprintf('Primaries + Multiples below surface %d',nl-1); 
   title(mytext) 
   xlabel('offset')
   ylabel('time')

end

PM=P0;clear P0;	
figure
for nl=1:nlayers
   
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
   %figure(5+nl),
   %subplot(223),
   %ww=fftshift(w);kxx=fftshift(kx);ww(1)=ww(2);kxx(1)=kxx(2);
   %wigb(fftshift(W1),1,kxx,ww);title('propagator in kx-w');
   %xlabel('kx')
   %ylabel('freq')

   %subplot(224),
   %wigb(w1,1,h,t);title('propagator in x-t');
   %xlabel('offset')
   %ylabel('time')

   W1=fft2(w1);
%  P=(W1.*(P+R1).*W1);

	r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft2(r1);
   
   if nl>1
		r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft2(r2);
	elseif nl==1
   	r2=zeros(NF,NH);r2(1,NH/2)=r0;R2=fft2(r2);
	end;

   PP=PM./(1+R2.*PM);
	PM=conj(W1).*PP-R1;
   p=ifft2(PP);%p(1:50,:)=0;PP=fft2(p);
   pm=ifft2(PM);%pm(1:50,:)=0;PM=fft2(pm);
   %p=ifft2(duplic2d(P));
   %figure(5+nl),
   
   subplot(221+2*(nl-1)),
   amx=max(max(abs(p)));
   wigb(p,1,h,t,amx);
   mytext=sprintf('Multiples from surface %d removed',nl-1); 
   title(mytext);
   xlabel('offset')
   ylabel('time')
   subplot(222+2*(nl-1))
   wigb(pm,1,h,t,amx);
   mytext=sprintf('Layer %d removed',nl); 
   title(mytext) 
   xlabel('offset')
   ylabel('time')

end
