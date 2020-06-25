function [LLO,LLU]=radonmat(w,h,p)
% Radon matrix LH*L
% such that d=Lm and m=LH*d;


nh=length(h);
np=length(p);
ii=1:np;
p=p(:).';
h=h(:);
F=exp(-i*w*((h.^2)*p));
FH=F';
LLO=FH*F;
LLU=F*FH;














