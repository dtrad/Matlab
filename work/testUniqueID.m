%test uniqueID
N=1000;
x = round((rand(1000,1)-0.5)*100);
y = round((rand(1000,1)-0.5)*100);

% for vectors with negative numbers use:
[index,rangey,minx,miny]=getUniqueId(x,y);
[x2,y2]=getXYfromId(index,rangey,minx,miny);

% for always positive numbers use:
%[index,rangey]=getUniqueId0(x,y);
%[x2,y2]=getXYfromId0(index,rangey);


more on
[index x x2 (x-x2) y y2 (y-y2)]
more off