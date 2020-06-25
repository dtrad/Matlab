function [FF,DD,VV]=mapC2R(F,d,v)
% Given the complex system of equations d=F*m
% returns the real matrices FF,DD, such that
% DD=FF*MM
% where MM=[real(m(:);imag(m(:))]
% Daniel Trad- September 1999
		% References: Numerical recipes, Press et al. 
FF=[real(F);imag(F)];                
FF2=[-imag(F);real(F)];
FF=[FF FF2];
DD=[real(d(:));imag(d(:))];
VV=[real(v(:));imag(v(:))];


