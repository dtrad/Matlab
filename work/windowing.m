function [W]=windowing(W)
% [W]=windowing(W)
% Daniel Trad- UBC.
[NF,NH]=size(W);
wind=hanning(NF);
wind=wind(:)*ones(1,NH);
W=W.*wind;


