function [u,p,h]=taup2xt_inv(v,hyperpar,sigma,option1,h,p,dt,Qp,Cn)
% 	[u,p,h]=taup2xt_inv(v,hyperpar,sigma,option1,h,p,dt,Qp,Cn)
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
% Daniel Trad-- 12-08-98
 
v=seis_shape(v);
[nt np]=size(v);
nh=np-5;

if nargin<7|isempty(dt) dt=0.004;end
if nargin<6|isempty(p) dp=1e-5;p=(0:np-1).*dp;end
if nargin<5|isempty(h) dh=25;h=(0:nh-1).*dh;end
if nargin<4|isempty(option1) option1='linear';end
if nargin<3|isempty(sigma) sigma=0.1;end
if nargin<2|isempty(hyperpar) hyperpar=10;end

np0=length(p);if np~=np0 display('offset axis does not match data');end

if option1(1:3)=='lin' Power=1;end
if option1(1:3)=='par' Power=2;dp=1e-8;p=(0:np-1).*dp;end
[message]=alias_check(p,h,60)
FS=1/dt;
w=2*pi*(0:nt-1).*FS/nt;

np=length(p);
nf=length(w);
nh=length(h);

p=p(:);
h=h(:).';

Cn=ones(np,1).*sigma;
Cn=diag(Cn);
CI=inv(Cn);

Qp=ones(nh,1)./hyperpar;
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

V=fft(v);
V=seis_shape(V);

for f=1:nt/2;
   % Operators
   F=exp(i*w(f)*(p*(h.^Power)));
   FH=F';
   L=F*WU;
   LH=FH*WV;
   if np > nh   
         u=inv(Qp+LH*CI*L)*LH*CI*(V(f,:)).';
   else    
         u=QpI*LH*inv(Cn+L*QpI*LH)*(V(f,:)).';
   end
   if(f==1) U=u;else U=[U,u];end,
end;
U=seis_shape(U);
U=duplic(U);
u=ifft(U);
