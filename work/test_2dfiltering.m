
close all
clear all

scale=1;
freq =0.01;
dx=25;
dy=50;
nx=64;
ny=32;
Nx=1./(2*dx);
Ny=1./(2*dy);

x=0:nx-1;x=x*dx;
y=0:ny-1;y=y*dy;y=y.';

%Assuming non separable axes generate the coordinates
X=ones(ny,1)*x(:).';
Y=y(:)*ones(1,nx);


Dx=sin(freq*2*pi*x);
Dy=sin(freq*1*pi*y);

Dx=Dx+sin(freq*1*pi*x);
Dy=Dy+sin(freq*0.5*pi*y);

d=Dy*Dx;
[ size(d) size(x) size(y) ]

filtery = 9:25;
filterx = 17:49;

FD = fft2(d);
FD2 = FD;
FD2(filtery,filterx)=0;

M = ones(size(FD2));
M2 = M;
M2(filtery,filterx)=0;

D  = ifft2(FD);
D2 = ifft2(FD2);

FD3 = fftshift(abs(fft2(D)));
FD4 = fftshift(abs(fft2(D2)));

close
figure(1)
subplot(221);
imagesc(d);

subplot(222)
plot(Dx);

subplot(223)
plot(Dy);

figure
subplot(221);
imagesc(abs(D));

subplot(222);
imagesc(abs(D2));

subplot(223);
imagesc(FD3);

subplot(223);
imagesc(FD4);

figure
mesh(FD3)

figure
mesh(FD4);

figure
imagesc(fftshift(M));

figure
imagesc(fftshift(M2))

figure
imagesc(FD3);

figure
imagesc(FD4);

