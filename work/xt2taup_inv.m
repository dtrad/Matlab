function [v,p,h]=xt2taup_inv(u,hyperpar,sigma,option1,h,p,dt)
% 		[v,p,h]=xt2taup_inv(u,hyperpar,sigma,option1,h,p,dt)
% Forward Transform through the inverse: 
% v=inv(Qp+L.CI.L*).L.CI.u
% u are the original data, 
% input 
%		u t-x data
% 		Qp regularization term
%     CI inverse of noise covariance matrix
% output 
%		v  tau-p data
%     % Optional:
%        h,  offset axis
%        p, p interval
%        dt, time interval
%        option1, linear or parabolic 
% Defaults:
% 		dt=0.004;
% 		dp=1e-5;p=(0:nh-1+10).*dp;
% 		dh=25;h=(0:nh-1).*dh;
% 		option1='linear';
% 		sigma=1e-3;
% 		hyperpar=1e-3;
%
% Daniel Trad-- 12-08-98

 
u=seis_shape(u);
[nt nh]=size(u);

if nargin<7|isempty(dt) dt=0.004;end
if nargin<6|isempty(p) dp=1e-5;p=(0:nh-1+10).*dp;end
if nargin<5|isempty(h) dh=25;h=(0:nh-1).*dh;end
if nargin<4|isempty(option1) option1='linear';end
if nargin<3|isempty(sigma) sigma=1e-3;end
if nargin<2|isempty(hyperpar) hyperpar=1e-3;end

nh0=length(h);if nh~=nh0 display('offset axis does not match data');end

if option1(1:3)=='lin' Power=1;end
if option1(1:3)=='par' Power=2;dp=1e-8;p=(0:nh-1+10).*dp;end



FS=1/dt;
w=2*pi*(0:nt-1).*FS/nt/2;

np=length(p);
nf=length(w);
nh=length(h);

p=p(:);
h=h(:).';

Cn=ones(nh,1).*sigma;
Cn=diag(Cn);
CI=inv(Cn);

Qp=ones(np,1).*hyperpar;
Qp=diag(Qp);
QpI=inv(Qp);

ii=1:np-1;
WV(ii)=p(ii+1)-p(ii);
WV(np)=WV(np-1);
WV=diag(WV);

ii=1:nh-1;
WU(ii)=h(ii+1)-h(ii);
WU(nh)=WU(nh-1);
WU=diag(WU);

UH=fft(u);
UH=seis_shape(UH);

for f=1:nt/2;
   % Operators
   F=exp(i*w(f)*(p*(h.^Power)));
   FH=F';
   L=F*WU;
   LH=FH*WV;
   if np < nh   
         v=inv(Qp+L*CI*LH)*L*CI*(UH(f,:)).';
   else    
      	QpI=diag(1./diag(Qp));
         v=QpI*L*inv(Cn+LH*QpI*L)*(UH(f,:)).';
   end
   if(f==1) V=v;else V=[V,v];end,
end;
V=V.';
V=duplic(V);
v=ifft(V);
