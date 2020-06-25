% Interferometry  
% Daniel Trad - CGGVeritas
clear
close all

% cube of data generated with fm_b_xw_m6.m
%load P11nomult; % no multiples
load P11;        % first order multiples
[NH,NH,NF]=size(P11);
NSG=32; % shot in the middle of the survey
NNUL=NF-100;
dt=0.004;t=(0:2*NF-1)*dt;
mute=0.3;
dh=15;hh=0:NH-1;hh=hh*dh;
f=70;wav=ricker(f,dt);wav=padzeros(wav,2*NF);

% bandpass data
S=fft(wav);
S=S(:);
S=[S(1:100);zeros(2*NF-200,1);S(2*NF-100+1:2*NF)];
SS=S(:)*ones(1,NH);

for ih=1:NH
    temp = P11(:,ih,:);
    temp=permute(temp,[3 1 2]);
    P11(:,ih,:)=(temp.*(SS(1:NF,:))).';    
end
    
[p1]=dsw2xt(P11(:,NSG,:));

imute = round(mute/dt);
MUTE=zeros(size(P11));
for ih=1:NH
    temp=dsw2xt(P11(:,ih,:));    
    temp(imute:end,:) = 0;
    temp = fft(temp);
    temp = temp.';
    
    MUTE(:,ih,:)=temp(:,1:NF);
end
    
P11 = P11 - MUTE;

OUT=zeros(size(P11));
NR = NH;
NS = NH;
for iw=1:NF-NNUL;
    for ir=1:NR
        for jr=1:NR
            temp = 0;
            for is=1:NS 
                temp = temp + P11(is,ir,iw)*conj(MUTE(is,jr,iw));
            end
            OUT(ir,jr,iw)=temp;
        end
    end
end;
P11=OUT;
clear OUT;
P11(:,:,NF-NNUL+1:NF)=0;
[p]=dsw2xt(P11(:,NSG,:));
p=ifft(fft(p).*SS);
%p=normalize(p)*coeff(1);

figure(1),
subplot(121);wigb(p1(1:2*NF,:),1,hh-hh(NSG),t);
title('Data shot 32(a)');ylabel('receiver');xlabel('time');

subplot(122),wigb(p(1:2*NF,:),1,hh-hh(NSG),t);
title('virtual shot 32(b)');xlabel('virtual receiver');ylabel('time');

%for ii=1:3:NH,figure,wigb(dsw2xt(P11(:,ii,:)));end
