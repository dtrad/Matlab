% Lab #2
% Exercise1 (c31.m)
% M is the number of coefficients to be computed

T0=pi;
N0=256;
T=T0/N0;
w0=2*pi/T0;

M=10;
t=0:T:T*(N0-1);t=t';

f=exp(-t/2);f(1)=0.604;
Dn=fft(f/N0)
[Dnangle,Dnmag]=cart2pol(real(Dn),imag(Dn));
k=0:length(Dn)-1;k=k';
figure(1)
subplot(211),stem(k,Dnmag);title('amplitude of FFT')
subplot(212), stem(k,Dnangle);title('Phase in Radians');

%Computation of trigonometric Fourier coefficients:

%c31;clg
C0=Dnmag(1); Cn=2*Dnmag(2:M);
Amplitudes=[C0;Cn]
theta0=Dnangle(1);theta=Dnangle(2:M);
Angles=Dnangle(1:M);
Angles=Angles*(180/pi);
disp('Amplitudes Angles')
[Amplitudes Angles]
k=0:length(Amplitudes)-1; k=k';
figure(2)
subplot(211),stem(k,Amplitudes);
title('Amplitude: First M Fourier coeffcients');
subplot(212), stem(k,Angles);
title('Phase: First M Fourier coeffcients');

%[Dnangle,Dnmag]=cart2pol(real(Dn),imag(Dn));
%Dn2=pol2cart(Dnangle,Dnmag);

f4=zeros(size(f));
f4=f4+C0;
figure,
subplot(5,2,1),plot(t,f4);
for ii=1:9 %length(Dnmag);
  f4=f4+Cn(ii)*cos(theta(ii)+ii*w0*t);
  subplot(5,2,ii+1),plot(t,f4);
end,



pause 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using duplic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
M=10
Dn2=[Dn(1:M);zeros(size(Dn(M+1:N0)))];
Dn3=duplic(Dn2(1:N0/2));
f4=N0*ifft(Dn3);

figure(3)
subplot(211),plot(t,real(f4));
title('f(t) reconstructed with the first 10 coefficients');
subplot(212), plot(t,f-real(f4));
title('f(t)-f4(t)');

M=10
Dn2=[Dn(1:M);zeros(size(Dn(M+1:N0)))];
Dn2(1)=0;
Dn3=duplic(Dn2(1:N0/2));
f4=N0*ifft(Dn3);

figure(4)
subplot(211),plot(t,real(f4));
title('f(t) reconstructed with the first 10 coefficients except the first one');
subplot(212), plot(t,f-real(f4));
title('f(t)-f4(t)');











