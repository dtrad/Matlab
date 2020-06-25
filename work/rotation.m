function [rot]=rotation(theta)
% function [rot]=rotation(theta)
% Function to rotate a vector 
% Seismology course eosc 354
% Daniel Trad - UBC 

rot(1,1)=cos(theta);
rot(1,2)=sin(theta);
rot(2,1)=-sin(theta);
rot(1,1)=cos(theta);