close all

load output
option1='parabolic';
hyperpar=1e-3;
sigma=1e-3;

xxc=xxc(1:512,1:70);figure,wigb(xxc)
[xxci,p,h]=xt2taup_inv(xxc,hyperpar,sigma,option1);figure;wigb(real(xxci));
xxcr=taup2xt(xxci,option1,h,p);figure,wigb(real(xxcr));
