function [theta]=max_angle(dt,dh,vel)
% [theta]=max_angle(dt,dh,vel)
% Computes max angle to avoid aliasing  for a given set of field parameters
% Daniel Trad- UBC.
theta=asin(dt*vel/dh)*180/pi;
