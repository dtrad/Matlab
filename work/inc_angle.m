function [theta]=inc_angle(dz,dh,NH)
% [theta]=inc_angle(dz,dh,NH)
% Angle of incidence for depth dz, offset h
% Daniel Trad- UBC.

h=(0:NH-1).*dh;
theta=acos(dz./((dz^2+h.^2).^0.5));
