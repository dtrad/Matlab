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
vel(2)=2000;
dt=0.004;
dh=25;
hh=0:NH-1;hh=hh*dh;
f=70;
nlayers=length(vel);
xa=dh*NH/2;
%xa=0;
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
%wav=[zeros(NF,NH/2) wav(:) zeros(NF,NH/2-1)];
S=fft(wav);
P=zeros(NF,NH);

for nl=nlayers:-1:2
   %[kz,kx,w]=obliquity(dt,dh,NF,NH,vel(nl));
	%kz=kz(1:NF,1:NH);
	%[kz,theta,k]=kz_fxz(dz,dh,NH,w,vel);
   r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft(r1);
   for ii=1:NH
        	[W1]=propagator(dt,dh,NF,NH,vel(nl),dz(nl),xa);
	end
	%W1=exp(-i*kz*dz(nl));
   w1=ifft(W1(:,:,1));
   
   %w1(1:NF/2,:)=0;
   %w1(:,NH/2+1:NH)=0;
   %W1=fft(w1);
   %w1=[w1 zeros(size(w1))];
   %w1=[w1;zeros(size(w1))];
   W12=fft2(w1);
   P2=W12.*fft2(r1);
   %.*W12;
   %W1=fft2(w1);
   for ii=1:NF;
   %P=(W1(1:NF,1:NH)*((P+R1).'));%.'*W1);
   P1(ii,:)=conv(W1(ii,:),R1(ii,:));
	end
   p1=ifft(P1);
   figure,
   subplot(211);
   wigb(p1)
   %p1(:,1:2*NH/4+30)=0;
   %p1(:,2*3*NH/4-30:2*NH)=0;
   subplot(212),
   wigb(p1)
	P1=fft(p1);
	for ii=1:NF
   %P11(ii,:)=P1(ii,:);
   P11(ii,:)=conv(W1(ii,:),P1(ii,2*NH/4+1:2*3*NH/4));
   end
end   
   
for ii=1:min(size(P11));
   P1(:,ii)=P11(:,ii).*S(:);
end
for ii=1:min(size(P2));
   P2(:,ii)=P2(:,ii).*S(:);
end
p2=ifft2(P2);

%P1=P11;
%P=P./(1+0*R1.*P);
%if nl==1 P=P.*S;end

p1=ifft(P1);
figure,
subplot(221),wigb(real(w1));
subplot(222),wigb(real(r1));
subplot(223),wigb(real(p1));
subplot(224),wigb(real(p2));

   



