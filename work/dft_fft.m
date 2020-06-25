function [D,kx,ky]=dft_fft(d,x,y,regx,regy);
% [D,kx,ky]=dft2(d,x,y);

if (nargin < 5) regy=0;end
if (nargin < 4) regx=0;end

dimx = 'irregul';
dimy = 'irregul';

if (regx==1) dimx='regular';end
if (regy==1) dimy='regular';end    

nx=length(x);
ny=length(y);
[ndy, ndx] = size(d);
if (ndx ~= nx)|(ndy ~= ny)
  display('error, size of d, x and y are not consistent');
  return;
end

  
dx=(x(end)-x(1))/(length(x)-1);
dy=(y(end)-y(1))/(length(y)-1);

kx=(-nx+2)/2:nx/2;
vv=fftshift(kx);kx=[vv(end) vv(1:end-1)]
kx=kx/(dx*nx);
kx=kx*2*pi;

ky=(-ny+2)/2:ny/2;
vv=fftshift(ky);ky=[vv(end) vv(1:end-1)];
ky=ky/(dy*ny);
ky=ky*2*pi;


%w=-lw/2:lw/2-1;
%w=w/(lw-1)*2*pi;
D=zeros(size(d));
if ( dimx == 'regular')
    for iy=1:ny
        D(iy,:)=fft(d(iy,:));
    end
else
    display('using dft along x');
    F=exp(i*(kx(:)*x(:).'));
    for iy=1:ny
        mv=d(iy,:);mv=mv(:);
        dv=F*mv;
        D(iy,:)=dv.';
    end
end

if ( dimy == 'regular')
    for ix=1:nx
        D(:,ix)=fft(D(:,ix));
    end
else
    display('using dft along y');
    F=exp(i*(ky(:)*y(:).'));
    for ix=1:nx
        D(:,ix)=F*D(:,ix);
    end
end

%keyboard

return
