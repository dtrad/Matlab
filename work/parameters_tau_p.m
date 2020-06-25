function [HH,WU0,WU1]=parameters_tau_p(h0,h1)
% This function receives a complete data set and returns
% a reduced data set to try interpolation and extrapolation.
% Given the original offset h0 and the final required offset h1
% computes the parameters HH, WU0, WU1, requiered by forwardnl3 and 
% the shrinked data u
%
% 	[u,HH,WU0,WU1]=seltrace(x,h0,h1)
%
% Input
% 		x original data  (REQUIRED)
%     h0 offset axis for forward transform
%     h1 offset axis for backward transform
%     
% output 
%   	u are the data to use
%     HH vector with 1 for valid traces 0 for missed traces.
%     WU0 diagonal matrix with offset interval for input data
%     WU1 diagonal matrix with offset interval for output data
%
% Daniel Trad-- UBC. Canada. 2-07-98

if nargin < 2 h1=h0;end;

h0=h0(:).';
h1=h1(:).';
nh0=length(h0);
nh1=length(h1);

HH=zeros(size(h1));

for ii=1:nh0;
	ind=find(h1==h0(ii));
	HH(ind)=1;
end

ii=1:nh0-1;
WU0=h0(ii+1)-h0(ii);
WU0(nh0)=WU0(nh0-1);
WU0=diag(WU0);

ii=1:nh1-1;
WU1=h1(ii+1)-h1(ii);
WU1(nh1)=WU0(nh1-1);
WU1=diag(WU1);




