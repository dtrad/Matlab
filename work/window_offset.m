function [W]=windo_offset(W)
% [W]=window_offset(W)
% Applies a hanning window along the offset traces
% Daniel Trad- UBC.
[NF,NH]=size(W);
wind=hanning(NH/2);wind=wind(:).';
wind=[wind(1:NH/4) ones(1,NH/2) wind(NH/4+1:NH/2)];
wind=ones(NF,1)*wind(:).';
W=W.*wind;


