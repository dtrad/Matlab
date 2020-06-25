% GEOP 520B-- Prof.: Tad Ulrych
% 	      Student: Daniel Trad
% Program to produce a band frequency wavelet.
clear

MM=128;
t=0:MM-1;
ts=-MM/2:MM/2-1;
f=(-MM/2:MM/2-1).*(1/MM);
w=zeros(size(t));
w(1)=1;
Y=fft(w);

figure(1);
subplot(221),plot(ts,fftshift(real(w)));title('Real spike- full band frequency');
subplot(222),plot(ts,fftshift(imag(w)));title('Imag spike- full band frequency');
subplot(223),plot(ts,fftshift(real(Y)));title('Real FT spike- full band frequency');
subplot(224),plot(ts,fftshift(imag(Y)));title('Imag FT spike- full band frequency');
figure(1)

figure(2);
subplot(211),semilogy(f,fftshift(abs(Y)));title('Spectrum- half band');
subplot(212),plot(f,fftshift(angle(Y)));title('Phase spectrum- half band');
figure(2);

YZ=zeros(size(Y(1:MM)));
%wind=hanning(MM/4);
wind=window('kais',MM/4);
IN=MM/8;INN=MM/4;
YZ(1)=0;
YZ(1+IN:IN+INN)=wind';
YZ(MM/2+1)=0;
YZ(MM/2+IN+2:MM/2+IN+INN+1)=wind';
%YZ=Y.*zwind;

wz=ifft(YZ(1:MM));

figure(3)
subplot(211),plot(ts,fftshift(real(wz)));title('Real spike constructed with Kaiser');
subplot(212),plot(ts,fftshift(imag(wz)));title('Imag spike constructed with Kaiser');
figure(3)
temp=fftshift(real(wz));
save c:\daniel\synth\kaiserspike temp

figure(4)
subplot(221),plot(f,fftshift(real(YZ(1:MM))));title('Real FT spike');
subplot(222),plot(f,fftshift(imag(YZ(1:MM))));title('Imag FT spike');
subplot(223),semilogy(f,fftshift(abs(YZ(1:MM))));title('mag spectrum');
subplot(224),plot(f,fftshift(angle(YZ(1:MM))));title('phase spectrum');
%gtext('spike constructed with KAISER window')
figure(4)

	

