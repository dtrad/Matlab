function [u,h0,h1,HH,WU0,WU1]=seltrace(xxc,nh0,nh1,dh,h_near,ntgap,bb,NZ)
% Select the traces to use.
%
% 	[u,h0,h1,HH,WU0,WU1]=seltrace(xxc,nh0,nh1,dh,h_near,ntgap,bb,NZ)
%
% This function receives a complete data set and returns
% a reduced data set to try interpolation and extrapolation.
% Input
% 		xxc original data 
% 		dh interval offset
%     h_near Closest trace in the original data
% 		nh0 Number of initial traces 
%		nh1 Number of final traces (nh0+n_interp)
% 		ntgap number of traces in the gap (traces missed)
% 		bb last trace before the gap
%     NZ Number of traces not used before the first
%     
% output 
%   	u are the data to use
%     h0 offset axis for forward transform
%     h1 offset axis for backward transform
%     HH vector with 1 for valid traces 0 for missed traces.
%     WU0 diagonal matrix with offset interval for input data
%     WU1 diagonal matrix with offset interval for output data
%
% Daniel Trad-- UBC. Canada. 2-07-98

if nargin < 7 NZ=0;end;
if nargin < 6 ntgap=0;end;
if nargin < 5 h_near=0;end;
if nargin < 4 dh=50;end;

% Index for u;

aa=1;

HH=ones(1,nh1);
HH(bb+1:bb+ntgap)=zeros(1,ntgap);
if nh0<nh1 HH(nh0+1:nh1)=zeros(1,nh1-nh0);end

u=shrinkt2(xxc(:,NZ+aa:NZ+nh0),HH);
h_near=h_near+NZ*dh;
clear temp;
nh=max(size(HH));
kk=0;
for ii=1:nh
   offset=h_near+(ii-1)*dh;
   if HH(ii)~=0 
      kk=kk+1;
      h0(kk)=offset;
      if kk>1 temp(kk-1)=h0(kk)-h0(kk-1);end
   end
   h1(ii)=offset;
end;   
temp(kk)=dh;

nh0=size(h0);
WU0=diag(temp);
WU1=dh*eye(nh1);

