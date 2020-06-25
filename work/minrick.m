% Program minrick.m
% Generate a minimum phase ricker wavelet
% Daniel Trad UBC-30-07-98


clear;
myclean
MM=128;
dt=0.1;
df=1/dt;
t=0:dt:(MM-1)*dt;
t0=6.4;
freq=(-MM/2:(MM/2-1)).*(df/MM);
w=6.2832;
x=(1-(w^2)*((t-t0).^2)/2).*exp(-(w^2)*((t-t0).^2)/4);

figure
plot(t,x);title('Ricker wavelet')
xlabel('time'),ylabel('amplitude');

X=fft(x);B=abs(X);

figure;
subplot(211),loglog(freq,B);title('Ricker wavelet-spectrum'),
xlabel('freq'),ylabel('amplitude');
subplot(212),semilogx(freq,angle(X));title('Ricker wavelet-spectrum'), 
xlabel('freq'),ylabel('phase');
lb=length(B);
PHI=imag(hilbert(log(B))); % Hilbert Transform of ln(|X|)
temp=exp(log(B)-i*PHI);temp=temp(:);
temp=duplic(temp(1:lb/2)); % FT of real time series
xmp=ifft(temp);
figure,
subplot(221);plot(real(xmp));title('MP ricker wavelet');xlabel('time');ylabel('amplitude');
subplot(222);plot(imag(xmp));title('MP ricker wavelet');xlabel('time');ylabel('amplitude');
subplot(223),loglog(freq,abs(temp));title('Ricker wavelet-spectrum'),
xlabel('freq'),ylabel('amplitude');
subplot(224),semilogx(freq,angle(temp));title('Ricker wavelet-spectrum'), 
xlabel('freq'),ylabel('phase');


B=log(abs(fft(x)));lb=length(B);
B(lb/2+2:lb)=0;
B(2:lb/2)=2.*B(2:lb/2);
xmp2=ifft(exp(B));

figure,
subplot(221);plot(real(xmp2));title('MP ricker wavelet');xlabel('time');ylabel('amplitude');
subplot(222);plot(imag(xmp2));title('MP ricker wavelet');xlabel('time');ylabel('amplitude');
subplot(223),loglog(fftshift(exp(B)));title('Ricker wavelet-spectrum'),
xlabel('freq'),ylabel('amplitude');
subplot(224),semilogx(fftshift(angle(exp(B))));title('Ricker wavelet-spectrum'), 
xlabel('freq'),ylabel('phase');
