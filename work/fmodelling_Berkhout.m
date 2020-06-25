% Forward Model- (Berkhout- 1982)
% Applies Berkhout formulation to create a full 2D survey
% using propagators, starting from top of halfspace going up,
% adding one layer at a time, and creating the multiples by 
% scattering of the wavefield at the different interfaces.
% Daniel Trad - (1997) UBC
%
% This like Lippman-Schwinger equation for scattering:
% P = P0 + P0 V P
% P = P0 + P0 V (P0 + P0 V P)
% P = P0 + P0 V P0 + P0 V P0 V P
% and so on.
% P0 can be taken as the Green Function (background wavefield for unitary
% source).
% The interpretation is that the wavefield P incident in V is a source
% which multiplied by the Green function produces the wavefield which
% adds to the existing background wavefield.
% The Green function is calculatted as the action of a propagator (W V W) 
% going first down, reflecting and then going up again.
% In Berkhout formulation the propagator can be expressed in (w,k) directly
% making it easier to convolve (just product in (w,k)) but to allow for
% horizontal variations it is better to express it in (w,x) and perform a 
% in space. In the code below the convolution is done through a matrix
% multiplication, putting the propagator in a circulant matrix.


clear
close all
NF=512;
NH=128;

minfreq = 1;
maxfreq = 125;
coeff(1)=1;
coeff(2)=1;
coeff(3)=1;

dz(1)=250;
dz(2)=350;
dz(3)=470;

vel(1)=2000;
vel(2)=2500;
vel(3)=3500;

%NULF=100; % Number of high frequencies which are not computed
r0=-1;
optionmult='first'
optionmult='none '
window=hanning(NH);
window=window*window.';
window=window+1e-2;

dt=0.004;t=(0:NF-1)*dt;
dh=15;
hh=0:NH-1;hh=hh*dh;
f=70;
Nyq = 1/(2*dt);
minIndex = round(minfreq/Nyq*(NF/2));
maxIndex = round(maxfreq/Nyq*(NF/2));

nlayers=length(vel);
xa=dh*NH/2; % shot in the center.
wav=ricker(f,dt);
wav=padzeros(wav,NF);
S=fft(wav);
P1=zeros(NF,NH);
P11=zeros(NH,NH,NF/2);
layers=coeff(:)*ones(1,NH);

%some tests
%layers(nlayers,1:NH/2-1)=0;
%layers(nlayers,NH/2+1:end)=0;
%layers(:,:)=0;
%layers(end,1)=1;

imagesc(layers);

% go layer by layer, freq by freq, starting from bottom.
% not sure why it uses the propagator at xa only..
for nl=nlayers:-1:1
   R1=diag(layers(nl,:));
   if (nl >1) R2=-diag(layers(nl-1,:));
   else R2=eye(NH)*r0;
   end
 
 
   % Create propagator (Green Function) from (xa,0) to (x,dz)
   [W1]=propagator(dt,dh,NF,NH,vel(nl),dz(nl),xa);
   w1=ifft(W1);w1=w1.*(ones(NF,1)*hanning(NH).');
   w1=normalize(w1)*0.13;
   W1=fft(w1);
   
   for ii=minIndex:maxIndex;
        %*********************************************************
        % Create a circulant matrix to perform the convolution of the 
        % propagator with the wavefield and reflectivity.        
        W11=[];
        W1=W1(:,NH:-1:1); % not sure if this flip is necessary
        for jj=1:NH/2;
            W11=[W11;W1(ii,NH/2-jj+1:NH) zeros(1,NH/2-jj) ];
        end
        for jj=1:NH/2;
            W11=[W11;zeros(1,jj) W1(ii,1:NH-jj)];
        end
        %*********************************************************
        % temporary wavefield (zero in first layer).
        temp=P11(:,:,ii); 
        % propagate wavefield from (xa,0) to (x,dz)
        % R1 acts as a source. 
        % This is like the scattering equation (P1 + P2 + ..)
        % where P1= W11*temp*W11 is existing wavefield going down and up.
        %       P2= W11*R1*W11 is the reflection in the interface at R1.
        temp=(W11*(temp+R1))*W11.';
        
        if optionmult=='inver'
            if max(max(abs(temp)))>1e-10
                if (cond(1-R2*temp+eye(NH)*0.05))<1e5
                    temp=temp*inv(1-R2*temp+eye(NH)*0.05);
                end
            end;
        elseif optionmult=='first'
            AMP=(1/2/(1e-3+max(max(abs(temp)))));
            temp=temp+AMP*R2*(temp.*temp);
            temp=temp+AMP^2*R2*R2*(temp.*temp.*temp);
        end;
        P11(:,:,ii)=temp;
   end
   P11(:,:,1:minIndex-1)=0;
   P11(:,:,maxIndex:NF/2)=0;   
   P1=P11(:,NH/2,:);
   P1=permute(P1,[3 1 2]);
   p1=ifft(duplic(windowing(P1(1:NF/2,:))));
   figure,wigb(p1);
   PF(:,:,nl)=P1;
end   
clear p   


for jj=1:nlayers   
	for ii=1:NH;
   	PF(:,ii,jj)=PF(:,ii,jj).*S(1:NF/2);
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

subplot(122),wigb(real(p(:,:,1)),1,hh-xa,t);title('2 layer and halfspace');xlabel('offset');ylabel('time')

figure,
wigb(real(p(:,:,1)),1,hh-xa,t);
title('Primaries- 2 layer model, (obtained with spatial convolution)')
xlabel('offset')
ylabel('time')

%P00=P11;
save P11.mat P11


