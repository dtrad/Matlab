function [c]=coeff(r1,r2,V1,V2)
% [c]=coeff(r1,r2,V1,V2)
% ri densities
% Vi velocities
% c reflection coefficients
c=(1/r1/V1-1/r2/V2)/(1/r1/V1+1/r2/V2);

