% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=64;
coeff(1)=0.2;
coeff(2)=0.3;
coeff(3)=0.5;
dz(1)=200;
dz(2)=250;
dz(3)=300;
vel(1)=1500;
vel(2)=2000;
vel(3)=2300;
NULF=0; % Number of high frequencies which are not computed
r0=-1;

dt=0.004;t=(0:NF-1)*dt;
dh=15;
hh=0:NH-1;hh=hh*dh;
f=70;
nlayers=length(vel)-1;
xa=dh*NH/2-1;
%xa=0;
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
%wav=[zeros(NF,NH/2) wav(:) zeros(NF,NH/2-1)];
S=fft(wav);
P1=zeros(NF,NH);
P11=zeros(NH,NH,NF/2-NULF);
for nl=nlayers:-1:1
   %r1=zeros(NF,NH);r1(1,NH/2)=coeff(nl);R1=fft(r1);
   R1=eye(NH)*coeff(nl);%R1(1,2:NH)=0;
   if nl>1
   R2=eye(NH)*coeff(nl-1);%R2(1,2:NH)=0;
   %r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
	elseif nl==1
   R2=eye(NH)*r0;
   %r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
	end

   
   [W1]=propagator(dt,dh,NF,NH/4,vel(nl),dz(nl),xa/4);
	w1=ifft(W1);
   W1=[zeros(NF,3*NH/4) W1];
   for ii=1:NF/2-NULF;
   %P=(W1(1:NF,1:NH)*((P+R1).'));%.'*W1);
   %P11(ii,:)=conv(W1(ii,:),(P1(ii,:)+R1(ii,:)));
   W11=[];
   W1=W1(:,NH:-1:1);
   for jj=1:NH;
      W11=[W11;zeros(1,jj-1) W1(ii,1:NH-jj+1)];
	end
%   for jj=1:NH/2;
%   W11=[W11;zeros(1,NH/2+jj-1) W1(ii,1:NH/2-jj+1)];
%	end
   
   temp=P11(:,:,ii);
	temp=(W11*(temp+R1))*W11.';
   %P11(ii,:)=(W11*P11(ii,:).').';%*W11.';
   %P11(ii,:)=P11(ii,:)*W11.';
	if max(max(abs(temp)))>1e-10
   %if (cond(1-R2*temp))<1e5 
    %  temp=temp*inv(1-R2*temp);
   %end
   end;
	P11(:,:,ii)=temp;

end
   
   P1=P11(:,NH/2,:);
   P1=permute(P1,[3 1 2]);
   P1=[P1;zeros(NULF,NH)];
   p1=ifft(duplic(windowing(P1(1:NF/2,:))));
   figure,wigb(p1);
   PF(:,:,nl)=P1;
end   
clear p   
%SS=S(:)*ones(1,NH);

for jj=1:nlayers   
	for ii=1:NH;
   	PF(:,ii,jj)=PF(:,ii,jj).*S(1:NF/2).';
	end
	p(:,:,jj)=ifft(duplic(windowing(PF(:,:,jj))));
end

figure,
subplot(221),wigb(real(w1),1,hh,t);title('propagator');xlabel('offset');ylabel('time')
t2=t(1:NF/2);
hh2=hh(NH/4+1:NH*3/4);

if nlayers > 2
   subplot(222),wigb(real(p(1:NF/2,NH/4+1:NH*3/4,3)),1,hh,t);title('1 layer and halfspace');xlabel('offset');ylabel('time')
end
if nlayers > 1
   subplot(223),wigb(real(p(1:NF/2,NH/4+1:NH*3/4,2)),1,hh2,t2);title('2 layer and halfspace');xlabel('offset');ylabel('time')
end

subplot(224),wigb(real(p(1:NF/2,NH/4+1:NH*3/4,1)),1,hh2,t2);title('3 layer and halfspace');xlabel('offset');ylabel('time')

figure,
wigb(real(p(1:NF/2,NH/4+1:NH*3/4,1)),1,hh2,t2(1:NF/2));
title('Primaries- 3 layer model, (obtained with spatial convolution)')
xlabel('offset')
ylabel('time')





