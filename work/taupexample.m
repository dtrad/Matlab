close all
load output
option1='parabolic';
xxc=xxc(1:512,1:40);figure,wigb(xxc)
[xxci,p,h]=xt2taup(xxc,option1);figure;wigb(real(xxci));
xxcr=taup2xt(xxci,option1,h,p);figure,wigb(real(xxcr));
