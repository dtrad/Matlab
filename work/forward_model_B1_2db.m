% Forward Model- (Berkhout- 1982)
% Daniel Trad - UBC
clear
close all
NF=256;
NH=32;
coeff(1)=0.5;
coeff(2)=0.5;
dz(1)=50;
dz(2)=150;
vel(1)=1000;
vel(2)=2000;
dt=0.004;
dh=10;
f=70;
nlayers=length(dz);
t=0:NF-1;t=t*dt;
h=0:NH-1;h=h*dh;
r0=-0.6;
%%%%%%%%%%%%%%%%%%%%%%
% Source
%
wav=rickerm(f,dt);
wav=padzeros(wav,NF);
for ii=1:NH
   wavm=[zeros(NF,ii-1) wav(:) zeros(NF,NH-ii)];
   %wavm=[wav(:) zeros(NF,NH-1)];
   S3(:,:,ii)=fft(wavm);
end
S3=permute(S3,[2,3,1]);


%%%%%%%%%%%%%%%%%%%%%%%

% Layers
%
P=zeros(NH,NH,NF);
for nl=nlayers:-1:1
   r1=ones(1,NH);r1=r1*coeff(nl);
   R1=diag(r1);
    
   if nl>1
      r2=ones(1,NH);r2=r2*coeff(nl-1);
      R2=diag(r2);
	elseif nl==1
      r2=ones(1,NH);r2=r2*r0;
      R2=diag(r2);
	end;
 

   [kz,kx,w,aa,alfa]=obliquity(dt,dh,NF,NH,vel(nl));
   tempW1=exp(-i*kz*dz(nl));
   %kk=(abs(kz)>1e-10);
   tempW1=tempW1.*aa;
   filter_alfa=((23*pi/180)>abs(real(alfa)));
   figure,wigb(filter_alfa);title('filter_alfa')
   tempW1=tempW1.*filter_alfa;
   %tempW1=zeros_2d(tempW1); % Clean the right dipping events
   w1=ifft2(tempW1);
   figure,wigb(w1);
   title('Propagator x,t')
   %w1=[w1 zeros(size(w1))];
   %w1(:,NH/2+1:NH)=0;
   %w1(NF/2+1:NF,:)=0;
   tempW1=fft(w1);
   figure,wigb(tempW1);title('Propagator in x,w')
   clear w1;
	for ii=1:NH
      W1(:,:,ii)=tempW1;
   end
   clear tempW1;
   W1=permute(W1,[2,3,1]);
   for ii=1:NF
      P(:,:,ii)=(W1(NH:-1:1,NH:-1:1,ii)*(P(:,:,ii)+R1)*W1(:,:,ii).');
      %P(:,:,ii)=P(:,:,ii).*(1-(P(:,:,ii).*R23(:,:,ii)));
      if nl==1 P(:,:,ii)=P(:,:,ii)*S3(:,:,ii);end
   end
   clear W1;
end
%%%%%%%%%%%%%%%%%%%%%%

P=permute(P,[3,1,2]);
for ii=1:NH
   p(:,:,ii)=ifft(P(:,:,ii));
   %p(:,:,ii)=ifft2(duplic2d(P(:,:,ii)));
end
figure
for ii=1:4:NH
   xx=round(ii/4);
   if xx<=3
      subplot(221+xx);
   elseif xx==4
      figure
      subplot(221);
   elseif xx>4   
      subplot(221+xx-4);
   end
   wigb(p(:,:,ii),1,h-(ii-1)*dh,t);
   title('synthetic data')
   xlabel('offset')
   ylabel('time')
   
end
