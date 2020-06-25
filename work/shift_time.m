function [y]=shift_time(x,threshold)
% [y]=shift_time(x,threshold)
% Given x, shift the time origin at the first x_i greater than thershold
% Default threshold=max(abs(x))./10;
% Daniel Trad-- UBC- 7-08-98

if nargin<2 threshold=max(abs(x))./10;end

lx=length(x);
ind1=min(find(abs(x)>threshold));
y=x(ind1:lx);
