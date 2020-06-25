clear
load c:\albert\zz;
close all
load output;
[u,haxis,ttaxis,vr]=interp_taup2(z,[],[],[],32,32,10,10,0,1,'linear ',1000,60,10,1e-3);
