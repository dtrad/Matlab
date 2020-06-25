function [u]=seis_shape(u)
% [u]=seis_shape(u)  output is time axis as rows
[mm nn]=size(u);if mm<nn u=u.';end
