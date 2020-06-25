% Example with Tau_p interpolation
clear
close all;
load c:\albert\zz;
z=seis_shape(z);
p=-0.004:1e-4:0.004;
[v,p,h]=xt2taup_inv(z,1,1e-3,'linear',[],p);
ur=taup2xt(v,'linear',h,p,0.004);
figure,
subplot(221),wigb(z);
subplot(222),wigb(v);
subplot(223),wigb(ur);
