function [D,dr,kx,ky]=dft2_irreg(d,x,y,dx,dy,nkx,nky,nxr,nyr,xr,yr,iter);
% [D,kx,ky]=dft2(d,x,y);
% This dft assumess input totally irregular so there is no differentiation
% along two axes. The operator has the total size (nx+ny) * (nx+ny) and it is applied
% once, rather than two smaller operators each one
% all the inputs are shaped as long vectors.

d=d(:);

nx=length(x);
ny=length(y);
nd = length(d);

if (nx ~= ny)|(nx ~= nd)
  display('error, size of d, x and y are not consistent');
  return;
end

kx=(-nkx+2)/2:(nkx-0)/2;
kx=kx/(dx*nkx);
kx=kx*2*pi;

ky=(-nky+2)/2:(nky-0)/2;
ky=ky/(dy*nky);
ky=ky*2*pi;

%w=-lw/2:lw/2-1;
%w=w/(lw-1)*2*pi;

%keyboard;

%generate unique k axis from the two kx, ky axes
KX=kx(:)*ones(1,nky);
KY=ones(nkx,1)*ky(:).';

K=[KX(:) KY(:)];

F=exp((i*(K(:,1)*x(:).'))+(i*(K(:,2)*y(:).')));
nn=length(KX(:));

D = F*d;
D = wtcgls(F',diag(1./D+1e-3),d,iter,0,1);

clear F;

XR=xr(:)*ones(1,nyr);
YR=ones(nxr,1)*yr(:).';
X=[XR(:).'; YR(:).'];
F=exp((i*(K(:,1)*X(1,:))+(i*(K(:,2)*X(2,:)))));

dr=real(F'*D);

% reconstruct to the whole size
% dr = ifft2(fftshift(reshape(D.',nkx,nky)),nxr,nyr,'symmetric');

return
