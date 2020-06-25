function [u]=taup2xt(v,option1,h,p,dt)
% Backward Transform: u~=L*v or u~=F*WV.v
% 	[u]=taup2xt(v,option1,h,p,dt)
% u are the original data, 
% Necesary arguments
% input:
%			v tau-p data
% output 
% 			u x-t data 
% Optional:
%        h,  offset axis
%        p, p interval
%        dt, time interval
%        option1, linear or parabolic 
% Daniel Trad-- 12-08-98
 

v=seis_shape(v);
[nt np]=size(v);

if nargin<5|isempty(dt) dt=0.004;end
if nargin<4|isempty(p) dp=1e-5;p=(0:np-1).*dp;end
if nargin<3|isempty(h) dh=25;h=(0:np-1+10).*dh;end
if nargin<2|isempty(option1) option1='linear';end

np0=length(p);if np~=np0 display('p axis does not match data');end

if option1(1:3)=='lin' Power=1;end
if option1(1:3)=='par' Power=2;end

%[message]=alias_check(p,h,[])
FS=1/dt;
w=2*pi*(0:nt-1).*FS/nt;

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

V=fft(v);
V=seis_shape(V);

for f=1:nt/2;
   F=exp(i*w(f)*(p*(h.^Power)));
   FH=F';
   ut=(FH*WV)*(V(f,:)).';
   if(f==1) UT=ut;
   else UT=[UT ut];
   end
end;
UT=UT.';
UT=duplic(UT);
u=ifft(UT);
