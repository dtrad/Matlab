clear,
close all
load c:\daniel\thesis\output

dt=0.004;
dh=25;
freq=70;
lambda=1e-3;
NF=512;
NH=64;
wav=ricker(freq,dt);
Dt=zeropadm(xxc(1:NF,1:NH),2*NF);

DTZ=zeros(size(Dt));
Dt=[Dt,DTZ];

% Redefine NF and NH because of the zero padding
[NF,NH]=size(Dt);

%figure,wigb(Dt);
%Dt=deconv_freq_domain_m(Dt,wav,lambda);

Dt(NF/2+1:NF,:)=0;  % Avoid wraparound
%figure,wigb(Dt);title('Data with wavelet deconvolved');

[D1t,D2t,kz]=born_multiples(Dt,dt,dh,wav);

D1t_old=D1t;
D2t_old=D2t;

D1t=deconv_freq_domain_m(D1t,wav,lambda);
D2t=deconv_freq_domain_m(D2t,wav,lambda);
D2t=deconv_freq_domain_m(D2t,wav,lambda);

scale=max(max(abs(Dt(150:512,:))))/max(max(abs(D1t(150:512,:))));
Dt_MR=real(Dt(1:512,1:64))+scale*real(D1t(1:512,1:64)+scale*real(D2t(1:512,1:64)));
figure,wigb(Dt_MR);title('Data-Multiple series')
figure,wigb(Dt(1:512,1:64));title('Data')
figure,wigb(D1t(1:512,1:64));title('Multiple series D1t')
figure,wigb(D2t(1:512,1:64));title('Multiple series D2t')
figure,wigb(D1t(1:512,1:64)+D2t(1:512,1:64));title('Multiple series D1t+D2t')