function [v]=forw_taup(u,h,p,dt,option1)
% Forward Transform: v=Lu or v=FWU.u
% u are the original data, 
% Necesary arguments
% input:
%			u t-x data
% output 
% 			V  tau-p data 
% Optional:
%        h,  offset axis
%        p, p interval
%        dt, time interval
%        option1, linear or parabolic 
% Daniel Trad-- 12-08-98
% 
u=seis_shape(u);
[nt nh]=size(u);

if nargin<5 option1='linear';end
if nargin<4 dt=0.004;end
if nargin<3 dp=1e-5;p=(0:nh+10).*dp;end
if nargin<2 dh=25;h=(0:nh-1).*dh;end

nh0=length(h);if nh~=nh0 display('offset axis does not match data');end

if option1(1:3)=='lin' Power=1;end
if option1(1:3)=='par' Power=2;end


FS=1/dt;
w=2*pi*(0:nt-1).*FS/nt/2;

np=length(p);
nf=length(w);
nh=length(h);

p=p(:);
h=h(:).';

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
   F=exp(i*w(f)*(p*(h.^Power)));
   v=(F*WU)*(UH(f,:)).';
   if(f==1) V=v;
   else V=[V v];
   end,   
end;
V=V.';
V=duplic(V);
v=ifft(V);
pause