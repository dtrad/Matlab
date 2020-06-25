%function [y,f]=shaping2d(epsilon,lf)
%Shaping filter
% epsilon regularization
% lf length of the filter to be used
% y: freq domain filtering
% y2: time domain filtering

close all
clear all

epsilon=1;

lx=128;ly=128;
ii=1:lx;
x=zeros(lx,ly);

x(lx/2-10,ly/2)=1;
x(lx/2+10,ly/2)=1;

w=zeros(size(x));

w(64,50:80)=1;
w(65,50:80)=1;

for i=1:20
  w(64+i-10,64+i-10)=1;
end

%w2=ricker(25,0.004);
w2=[0.1 -0.3 1 -0.3 0.2];
w=conv2(w,w2);
w=w(1:lx,1:ly);

simage(w);colorbar;

y=real(fftshift(ifft2(fft2(x).*fft2(w))));
%y=conv2(x,w);
simage(y);

X=fft2(x);
W=fft2(w);




%f=real(ifft(conj(fft(w)).*fft(x)./(epsilon+conj(fft(w)).*fft(w))));
%f=real(ifft(fft(x)./(epsilon+abs(fft(w)))));
F=(conj(W).*X./(max(epsilon,(conj(W).*W))));

%figure,simage(real(conj(W).*(W)))


Y2=W.*F;



Y2=duplic2d(Y2);
y2=real(ifft2(Y2));

F=duplic2d(F);
f=real(fftshift(ifft2(F)));

figure
subplot(221);simage(x);colorbar;title('signal before convolution')
%subplot(222);simage(w);colorbar;title('2D wavelet')
subplot(222);simage(y);colorbar;title('convolution y=x*w')
subplot(223);simage(y2);colorbar;title('2D deconvolution x=y*^{-1}w')
subplot(224);simage(f);colorbar;title('filter')

% Time domain 2d deconvolution
y3=conv2(y,f);

subplot(223);simage(y2);colorbar;title(['2D deconvolution x=y*^{-1}w')

figure 

return;

y2=conv(w,f(1:lf));
y2=y2(1:lx);

figure
subplot(311);plot(ii,x,ii,w,ii,f)
subplot(312);plot(ii,y)
subplot(313);plot(ii,y2)


