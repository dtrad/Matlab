% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=32;
coeff(1)=0.2;
coeff(2)=0.3;
coeff(3)=-0.5;
dz(1)=300;
dz(2)=550;
dz(3)=300;
vel(1)=1000;
vel(2)=2000;
vel(3)=2300;
r0=-0.8;

dt=0.004;
dh=15;
hh=0:NH-1;hh=hh*dh;
f=70;
nlayers=length(vel);
xa=dh*NH/2;
%xa=0;
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
%wav=[zeros(NF,NH/2) wav(:) zeros(NF,NH/2-1)];
S=fft(wav);
P1=zeros(NF,NH);
P11=zeros(NF,NH);
for nl=nlayers:-1:1
   %r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft(r1);
   R1=ones(1,NH)*coeff(nl);R1(1,2:NH)=0;
   if nl>1
   R2=ones(1,NH)*coeff(nl-1);R2(1,2:NH)=0;
   %r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
elseif nl==1
   R2=ones(1,NH)*r0;
   %r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
	end

   for ii=1:NH
        	[W1]=propagator(dt,dh,NF,NH,vel(nl),dz(nl),xa);
	end
   w1=ifft(W1);
   
   for ii=1:NF;
   %P=(W1(1:NF,1:NH)*((P+R1).'));%.'*W1);
   %P11(ii,:)=conv(W1(ii,:),(P1(ii,:)+R1(ii,:)));
   W11=[];
   for jj=1:NH;
   W11=[W11;W1(ii,jj:NH) W1(ii,1:jj-1)];
	end
   P11(ii,:)=(W11*(P1(ii,:).'+R1.')).';
   %P11(ii,:)=(W11*P11(ii,:).').';%*W11.';
   %P11(ii,:)=P11(ii,:)*W11.';
   %P11(ii,:)=P11(ii,:)*inv(1-R2(:)*P11(ii,:));
   P11T(ii,:)=P11(ii,:)*inv(1-ones(NH,1)*P11(ii,:));

end
P11=P11T;
   NHH=min(size(P11));
   
   %P1=P11(:,(NHH+1)/2-NH/2+1:(NHH+1)/2+NH/2);
   P1=P11;
   p1=ifft(P1);
   %p1=normalize(p1);p1=p1*coeff(nl);
   figure,wigb(p1);
   %P1=fft(p1);
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
subplot(221),wigb(real(w1));
subplot(222),wigb(real(p(:,:,1)));
subplot(223),wigb(real(p(:,:,2)));
subplot(224),wigb(real(p(:,:,3)));

figure,
wigb(real(p(1:NF/2,:,1)));



