function [D,kx,ky]=dft2b(d,x,y);
% [D,kx,ky]=dft2(d,x,y);

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

F=exp(i*(kx(:)*x(:).'));
size(F)
size(d)


D=zeros(size(d));



for ix=1:ny
  mv=d(ix,:);mv=mv(:);
  dv=F*mv;
  D(ix,:)=dv.';
end

F=exp(i*(ky(:)*y(:).'));

for iy=1:nx
  D(:,iy)=F*D(:,iy);
end

%keyboard

return
