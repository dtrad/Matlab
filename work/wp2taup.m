function [Rt]=wp2taup(R)
R=duplic(R);
Rt=ifft(R);
figure,wigb(real(Rt));
