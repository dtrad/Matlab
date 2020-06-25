function [Q,D]=attenuation(tau,p,V)
if nargin <3 V=2000;end
tau=tau(:);
TAU=tau*ones(1,length(p));
p=p(:).';
P=ones(length(tau),1)*p;
D=distance_taup(TAU,P,V);
Q=1./(1+D);
%Q=(1-(D./max(max(D))));
