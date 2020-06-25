% GEOP 527B
% Problem 6- Assignment 3.
% 11-11-97
% Daniel Trad

%	HILBERT(X) is the Hilbert transform of the real part
%	of vector X.  The real part of the result is the original
%	real data; the imaginary part is the actual Hilbert
%	transform.  
%	Author(s): C. Denham, 1-7-88
%		   L. Shure, 11-19-88, 5-22-90 - revised
%		   K. Creager, 3-19-92, modified to use power of 2 FFT
%		   T. Krauss, 11-4-92, revised
%	Reference(s):
%	  [1] Jon Claerbout, Introduction to Geophysical Data Analysis.

clear;
for i=1:20;close;end;
MM=128;
dt=0.1;
df=1/dt;
t=-(MM/2*dt):dt:(MM/2-1)*dt;
t=0:dt:(MM-1)*dt;
%t=0:dt:(MM-1)*dt;
t0=6.4;
freq=(-MM/2:(MM/2-1)).*(df/MM);
w=6.2832;
f=(1-(w^2)*((t-t0).^2)/2).*exp(-(w^2)*((t-t0).^2)/4);

R=1;
theta=30;

figure(1);
plot(t,f);title('Ricker wavelet'),
xlabel('time'),ylabel('amplitude');
figure(1);

F=fft(f);

figure(2);
subplot(211),plot(freq,fftshift(abs(F)));title('Ricker wavelet-spectrum'),
xlabel('freq'),ylabel('amplitude');
subplot(212),plot(freq,fftshift(angle(F)));title('Ricker wavelet-spectrum'), 
xlabel('freq'),ylabel('phase');
figure(2);

g=imag(hilbert(f));

for i=1:3;
figure;
 
theta1=theta*1+(i-1)*120;
theta2=theta*2+(i-1)*120;
theta3=theta*3+(i-1)*120;
theta4=theta*4+(i-1)*120;

omega1=theta1*pi/180;
omega2=theta2*pi/180;
omega3=theta3*pi/180;
omega4=theta4*pi/180;

rf1=R*(f.*cos(omega1)+g.*sin(omega1));
rf2=R*(f.*cos(omega2)+g.*sin(omega2));
rf3=R*(f.*cos(omega3)+g.*sin(omega3));
rf4=R*(f.*cos(omega4)+g.*sin(omega4));

text1=sprintf('theta=%d',theta1);
text2=sprintf('theta=%d',theta2);
text3=sprintf('theta=%d',theta3);
text4=sprintf('theta=%d',theta4);

subplot(221);plot(rf1);title(text1);xlabel('time');ylabel('amplitude');
subplot(222);plot(rf2);title(text2);xlabel('time');ylabel('amplitude');
subplot(223);plot(rf3);title(text3);xlabel('time');ylabel('amplitude');
subplot(224);plot(rf4);title(text4);xlabel('time');ylabel('amplitude');

end;

