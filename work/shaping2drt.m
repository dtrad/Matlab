function [f,F]=shaping2drt(d,dout)
%Shaping filter
% epsilon regularization
% lf length of the filter to be used
% y: freq domain filtering
% y2: time domain filtering

%close all
%clear all

%load sudata.mat

f1=50;
n1=128;
f2=1;
n2=64;

w=d(f1:f1+n1-1,f2:f2+n2-1);
x=dout(f1:f1+n1-1,f2:f2+n2-1);
y=w;

epsilon=1;
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
subplot(223);simage(y2);colorbar;title('2D deconvolution FD x=y*^{-1}w')
subplot(224);simage(f);colorbar;title('filter')

% Time domain 2d deconvolution
y3=conv2(y,f);

n12=n1/2;
n22=n2/2;

y3=y3(n12:n12+n1-1,n22:n22+n2-1);

figure
subplot(221);simage(x);colorbar;title('signal before convolution')
%subplot(222);simage(w);colorbar;title('2D wavelet')
subplot(222);simage(y);colorbar;title('convolution y=x*w')
subplot(223);simage(y3);colorbar;title('2D deconvolution TD x=y*^{-1}w')
subplot(224);simage(f);colorbar;title('filter')

% Apply to another event





return;



