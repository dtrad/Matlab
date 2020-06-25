function [d]=distance_taup(tau,p,V)
% [d]=distance_taup=(tau,p,V)
% Output d: distance travelled by a ray at tau and p;
% q=vertical slowness.
% After the critical angle I take the absolute value
% Daniel Trad - UBC - 
if nargin <3 V=2000;end
costheta=sqrt((1-(p*V).^2));
q=costheta/V;
d=tau./(q+eps);
