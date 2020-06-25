n=4000000;
n=100000;
%load ~/sugeomfix/p0.sg
%load ~/sugeomfix/p1.sg
load ~/modules/geomfix/p0.sg
load ~/modules/geomfix/p2.sg
p1=p2;
minx=375589
miny=6014537
minx=0
miny=0

figure(1),
plot(p0(1:n,1)-minx,p0(1:n,2)-miny,"+;original;",p1(1:n,1),p1(1:n,2),".;rotated;");
xlabel="shots"

figure(2),
plot(p0(1:n,3)-minx,p0(1:n,4)-miny,"+;original;",p1(1:n,3),p1(1:n,4),".;rotated;");
title="receivers"