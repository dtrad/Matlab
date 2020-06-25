function [del]=delay(delay0,p,v)
% function [del]=delay(delay0,p,v)
% Computes the delay as a function of p for the predictive deconvolution
% Input: 
%      delay0 : delay for p=0
%      p: p axis for every trace
%      v: velocity of the multiple.
% Output:
%      del(pp);

% Daniel Trad - Jun 4-98

np=max(size(p));

pp=1:np;
del(pp)=delay0*((1-(p(pp)*v).^2)).^2;




