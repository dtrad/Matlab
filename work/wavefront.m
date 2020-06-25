function [P]=wavefront(r1,r2,ks,kg,w)
xs=r1(1,:);
zs=r1(2,:);
xg=r2(1,:);
zg=r2(2,:);

P=exp(+i*(ks*(xg-xs)+kg*abs(zg-zs)));
