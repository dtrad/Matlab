function [r]=residual(alfa,d,B,m)
r=(d-alfa*B*m)'*(d-alfa*B*m);

