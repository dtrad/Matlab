function [y]=normalize(x)
% [y]=normalize(x) such that max(y)=1
% Daniel Trad- UBC
x=seis_shape(x);
y=x./max(max(abs(x)));
