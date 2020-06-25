function [zz]=triang(x,y,z,xx,yy)
% Example of point estimation by triangulation
% given three points with coordinates x,y,z
% find the value of (xx,yy,zz) using triangulation

x1=x(1);y1=y(1);
x2=x(2);y2=y(2);
x3=x(3);y3=y(3);

area1= area(xx,yy,x2,y2,x3,y3)
area2 = area(xx,yy,x1,y1,x3,y3)
area3 = area(xx,yy,x1,y1,x2,y2)

zz=z(1)*area1+z(2)*area2+z(3)*area3;
zz=zz/area(x1,y1,x2,y2,x3,y3);
return

function [d]=distance(x,y,xx,yy);
d=sqrt((x-xx)^2+(y-yy)^2);
return

function [a]=area(x1,y1,x2,y2,x3,y3)

d1=distance(x1,y1,x2,y2);
d2=distance(x1,y1,x3,y3);
a=0.5*sin(myangle(x1,y1,x2,y2,x3,y3,d1,d2))*d1*d2;

return;

function [ang]=myangle(x1,y1,x2,y2,x3,y3,d1,d2)
ang=acos(crossproduct(x1,y1,x2,y2,x3,y3)/(d1*d2));
return;

function [cp]=crossproduct(x1,y1,x2,y2,x3,y3)
cp=(x2-x1)*(x3-x1)+(y2-y1)*(y3-y1);
return;
