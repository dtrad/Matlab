% Test DFT in two dimensions for an irregularly sampled data set
% 1) Create an oversampled data set 
% 2) Decimate it in an irregular way
% 3) Compare DFT2 with FFT2 (the last one for the original data set)
clear all
close all
threshold = 0.25;
scale=1;
freq =0.01;
dx=25;
dy=50;
nx=64;
ny=32;
nkx=nx;
nky=ny;
iter=5;

Nx=1./(2*dx);
Ny=1./(2*dy);

x=0:nx-1;x=x*dx;x=x(:);
y=0:ny-1;y=y*dy;y=y(:);

example=2;
if (example == 1)
    Dx=sin(freq*2*pi*x);
    Dy=ones(1,ny);
    d=Dx*Dy;
elseif (example == 2)
    w=ricker(15,0.01);
    d=zeros(nx,ny);
    for ix=1:nx 
        iy=1*ix;    
        if (iy<ny)
            d(ix,iy)=1;
        end
    end
    for iy=1:ny;
        dw=conv(d(:,iy),w);
        d(:,iy)=dw(1:nx);
    end
end


[ size(d) size(x) size(y) ]
DD = fftshift(fft2(d));

[dd,xx,yy,nn,ds]=decimate_irreg(d,x,y,nx,ny,threshold);

[D,dr,kx,ky]=dft2_irreg(dd,xx,yy,dx,dy,nkx,nky,nx,ny,x,y,iter);
D = reshape(D,nkx,nky);
dr = reshape(dr,nx,ny);

close
figure(1)
subplot(221);
imagesc(d);

subplot(222)
%plot(Dx);

subplot(223);
imagesc(abs(D));

subplot(224);
imagesc(abs(DD));

figure(2)
mesh(abs(D))

figure(3)
mesh(abs(DD))

figure(4);imagesc(ds);

figure(5);imagesc(dr);
figure(6);imagesc(d);
figure(7);imagesc(d.*ds);




