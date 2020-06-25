function [h]=rand_offset(nh,dh0,scale)

% Computes an offset axis of inteval given by dh0+perturbation
% Daniel Trad- UBC- 22 March 1999
 
if nargin < 2 dh0=25;end
dhp=randn(nh,1);
h(1)=0;for i=2:nh;h(i)=scale*round((dhp(i)))+dh0+h(i-1);end;
%h=h-round((max(h))/2);
