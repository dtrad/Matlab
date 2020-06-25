n1=128;n2=64;
f1=0.4;
f2=0.15;
d1=1;
d2=1;
pert1 = 0.2;
pert2 = 0.5;
nfig=1;
% regular 
x1=0:n1-1;x1=x1*d1;
x2=0:n2-1;x2=x2*d2

% irregularity along x1
perturb1=1+pert1*(rand(n1,1)*2-1);
x1(1)=0;

for ix=2:n1
    x1(ix)=x1(ix-1)+d1*perturb1(ix-1);
end

% irregularity along x2
perturb2=1+pert2*(rand(n2,1)*2-1);
x2(1)=0;

for ix=2:n2
    x2(ix)=x2(ix-1)+d2*perturb2(ix-1);
end

figure(1); 
subplot(211);plot(1:n1,x1,'o-');
subplot(212);plot(1:n2,x2,'o-');

for i1=1:n1
    for i2=1:n2
        z(i1,i2)=sin(x1(i1)*2*pi*f1)+sin(x2(i2)*2*pi*f2);
    end
end


figure(2),
subplot(221);imagesc(z)

ZFFT=abs(fftshift(fft2(z)));
subplot(222);mesh(ZFFT);colorbar;

%ZDFT=abs(fftshift(dft2b(z,x2,x1)));
%figure,imagesc(ZDFT);colorbar;

ZDFT=abs(fftshift(dft_fft(z,x2,x1,0,0)));
subplot(223);mesh(ZDFT);colorbar;
    
subplot(224);imagesc(ZFFT-ZDFT);colorbar;

