scale=1;
freq =0.01;
dx=25;
dy=50;
nx=128;
ny=64;
Nx=1./(2*dx);
Ny=1./(2*dy);

x=0:nx-1;x=x*dx;
y=0:ny-1;y=y*dy;y=y.';

Dx=sin(freq*2*pi*x);
Dy=ones(ny,1);
d=Dy*Dx;
[ size(d) size(x) size(y) ]

[D,kx,ky]=dft2(d,x,y);

DD = fftshift(fft2(d));

close
figure(1)
subplot(221);
imagesc(d);

subplot(222)
plot(Dx);

subplot(223);
imagesc(abs(D));

subplot(224);
imagesc(abs(DD));

figure
mesh(abs(D))

figure
mesh(abs(DD))

