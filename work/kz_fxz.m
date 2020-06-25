function [kz,theta,k]=kz_fxz(dz,dh,NH,w,vel)
% [kz]=kz(dz,dh,NH,w,vel)
% Daniel Trad- UBC
k=w./vel;k=k(:);
theta=inc_angle(dz,dh,NH);
kz=k*cos(theta);