close all

load output
option1='linear';
hyperpar=1e-3;
sigma=1e-3;

xxc=xxc(1:512,1:40);figure,wigb(xxc)
[xxci,p,h]=xt2taup(xxc,option1);figure;wigb(real(xxci));
xxcr=taup2xt_inv(xxci,hyperpar,sigma,option1,h,p);figure,wigb(real(xxcr));
