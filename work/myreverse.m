function [u]=reverse(u)
% Reverse the time axis for a matrix
% Daniel Trad- UBC,
u=seis_shape(u);
lx=max(size(u));
u=u(lx:-1:1,:);
