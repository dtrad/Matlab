function [Gx,Gy,G]=grad2d(n,dx,dy)
if (nargin < 1) n=8;end
if (nargin < 2) dx=1;end
if (nargin < 3) dy=dx;end
% GRID for u

GRDn=1:1:(n+1)^2;
GRDn=reshape(GRDn,n+1,n+1);

% GRID for Jx

GRDex=1:1:(n+1)*n;
GRDex=reshape(GRDex,n,n+1);

% GRID for Jy

GRDey=1:1:(n+1)*n;
GRDey=reshape(GRDey,n+1,n);

ix=[];
jx=[];
kx=[];
iy=[];
jy=[];
ky=[];

ix=mkvc(GRDex);
jx=mkvc(GRDn(1:end-1,:));
kx=mkvc(-1/dx*ones(n*(n+1),1));

ix=[ix;ix];
jx=[jx;mkvc(GRDn(2:end,:))];
kx=[kx;-kx];

Gx=sparse(ix,jx,kx);

%Test
%spy(Gx)
%FGx=full(Gx);
%FGx(1:10,1:10)

iy=mkvc(GRDey);
jy=mkvc(GRDn(:,1:end-1));
ky=mkvc(-1/dy*ones(n*(n+1),1));

iy=[iy;iy];
jy=[jy;mkvc(GRDn(:,2:end))];
ky=[ky;-ky];

Gy=sparse(iy,jy,ky);

%spy(Gy)
%FGy=full(Gy);
%FGy(1:10,1:10)


G=[Gx;Gy];
