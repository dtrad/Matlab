% Krigging interpolation: example.
% c0: nugget effect
% c0+c1: sill
% a:  range;
load krigdata;
c0=0;
c1=10;
a=10;

minx=min(x);
maxx=max(x);
miny=min(y);
maxy=max(y);
index=0;
for j=miny:maxy
    for i=minx:maxx
        index=index+1;
        xx(index)=i;
        yy(index)=j;
    end
end
zz=kriging(c0,c1,a,x,y,z,xx,yy);
figure(1);plot3(x,y,z,'o');figure(gcf);
figure(2);plot3(x,y,z,'o',xx,yy,zz,'.');figure(gcf);
    
