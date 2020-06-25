function [D1t,D2t,kz]=born_multiples(Dt,dt,dh,wav)
% [D1t,t,h]=born_multiples(Dt,dt,dh,wav)
% Daniel Trad- UBC-CA

[NF,NH]=size(Dt);
D=fft2(Dt(1:NF,1:NH));
kz=obliquity(dt,dh,NF,NH,1000);

figure,wigb(fftshift(real(D)));title('real(FFT2(Data))')
figure,wigb(fftshift(imag(D)));title('imag(FFT2(Data))')

D1=D.*D*i/pi.*kz;

D1=duplic2d(D1);
D1t=ifft2(D1);
temp=D1t(1:NF/2,1+32:NH-32);


D1t=[temp;zeros(size(temp))];clear temp;
D1t=[D1t,zeros(size(D1t))];


D1=fft2(D1t);
D2=D.*D1.*kz*i/pi;

figure,
wigb(D1t)

display ('max imag in D1t=')
max(max(abs(imag(D1t))))

D2=duplic2d(D2);
D2t=ifft2(D2);
temp=D2t(1:NF/2,1+32:NH-32);
D2t=[temp;zeros(size(temp))];clear temp;
D2t=[D2t,zeros(size(D2t))];

figure,
wigb(D2t)

display ('max imag in D2t=')
max(max(abs(imag(D2t))))