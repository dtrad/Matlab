% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=512;
NH=64;
coeff(1)=0.3;
coeff(2)=0.3;
%coeff(3)=0.5;
dz(1)=200;
dz(2)=450;
%dz(3)=300;
vel(1)=2000;
vel(2)=3000;
%vel(3)=2300;
NULF=100; % Number of high frequencies which are not computed
r0=-0.8;
optionmult='nomul'
window=hanning(NH);
window=window*window.';
window=window+1e-2;

dt=0.004;t=(0:NF-1)*dt;
dh=15;
hh=0:NH-1;hh=hh*dh;
f=70;
nlayers=length(vel);
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
   	R2=eye(NH)*(-coeff(nl-1));%R2(1,2:NH)=0;
   	%r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
		elseif nl==1
   	R2=eye(NH)*r0;
   	%r2=zeros(NF,NH);r2(1,NH/2)=coeff(nl-1);R2=fft(r2);
	end

   
   [W1]=propagator(dt,dh,NF,NH,vel(nl),dz(nl),xa);
   w1=ifft(W1);w1=w1.*(ones(NF,1)*hanning(NH).');
   w1=normalize(w1)*0.13;
   %W1=fft(duplic(w1(1:NF/2,:)));
   W1=fft(w1);
   
   %W1=[zeros(NF,3*NH/4) W1];
   for ii=1:NF/2-NULF;
   	%P=(W1(1:NF,1:NH)*((P+R1).'));%.'*W1);
   	%P11(ii,:)=conv(W1(ii,:),(P1(ii,:)+R1(ii,:)));
   	W11=[];
   	W1=W1(:,NH:-1:1);
   	for jj=1:NH/2;
      	W11=[W11;W1(ii,NH/2-jj+1:NH) zeros(1,NH/2-jj) ];
		end
   	for jj=1:NH/2;
      	W11=[W11;zeros(1,jj) W1(ii,1:NH-jj)];
		end
   
%   for jj=1:NH/2;
%   W11=[W11;zeros(1,NH/2+jj-1) W1(ii,1:NH/2-jj+1)];
%	end
   
   	temp=P11(:,:,ii);
   %if nl==1;keyboard;end;
		temp=(W11*(temp+R1))*W11.';
   %P11(ii,:)=(W11*P11(ii,:).').';%*W11.';
   %P11(ii,:)=P11(ii,:)*W11.';
   	if optionmult=='inver'
			if max(max(abs(temp)))>1e-10
   			if (cond(1-R2*temp+eye(NH)*0.05))<1e5
      			temp=temp*inv(1-R2*temp+eye(NH)*0.05);
   			end
			end;
   	elseif optionmult=='first'
      	AMP=(8/NH/(1e-2+max(max(abs(temp)))));
      	%AMP=1;
      	%temp=window.*temp;
   
      	temp=temp+AMP*temp*R2*temp;
      	%temp=temp+AMP^2*temp*R2*temp*R2*temp;
      	%temp=temp+AMP.^3*temp*R2*temp*R2*temp*R2*temp;
      	%temp=temp./window;
   	end;
		P11(:,:,ii)=temp;

	end

	P11(:,:,NF/2-NULF:NF/2)=0;   
	P1=P11(:,NH/2,:);
	P1=permute(P1,[3 1 2]);
	%P1=[P1;zeros(NULF,NH)];
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

if nlayers > 2
   subplot(222),wigb(real(p(:,:,3)),1,hh-xa,t);title('1 layer and halfspace');xlabel('offset');ylabel('time')
end
if nlayers > 1
   subplot(223),wigb(real(p(:,:,2)),1,hh-xa,t);title('2 layer and halfspace');xlabel('offset');ylabel('time')
end

subplot(224),wigb(real(p(:,:,1)),1,hh-xa,t);title('3 layer and halfspace');xlabel('offset');ylabel('time')

figure,
wigb(real(p(:,:,1)),1,hh-xa,t);
title('Primaries- 3 layer model, (obtained with spatial convolution)')
xlabel('offset')
ylabel('time')

if optionmult=='first'
save c:\daniel\thesis\P11.mat P11
elseif optionmult=='nomul'
P00=P11;
save c:\daniel\thesis\P00.mat P00
end

