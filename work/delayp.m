function [del]=delayp(delay0,p,v,dt)
% function [del]=delay(delay0,p,v)
% Computes the delay as a function of p for the predictive deconvolution
% Input: 
%      delay0 : delay for p=0 in sec.
%      p: p axis for every trace
%      v: velocity of the multiple.
% Output:
%      del(pp); delay time in sec.

% Daniel Trad - Jun 4-98

if nargin<4 dt=0.004;end;

np=max(size(p));

pmax=max(find(p<1/v));
pp=1:pmax-1;

%delay0=(ddelay0-1)*dt;
del(pp)=delay0*((1-(p(pp)*v).^2)).^2;
del(pmax:np)=0;%zeros(pmax:np);
%ddel(1:np)=round(del(1:np)./dt)+1; delay in index



