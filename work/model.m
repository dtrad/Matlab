function[Bx,By,Bz,xspacing,yspacing]=model(numdipoles)

xspacing=input('Enter the spacing of x data points. \n');
yspacing=input('Enter the spacing of y data points. \n');
c=input('Enter the length of your x data array. \n');
d=input('Enter the length of your y data array. \n');

for i=1:numdipoles

   if i~=1
      Bxtemp=Bx;
      Bytemp=By;
      Bztemp=Bz;
   end

xs=input(sprintf('Enter the x coordinate for dipole #%g. \nKeep in mind that your survey arbitrarily starts from \ncoordinates (0,0) and works in the positive direction. \n',i));
ys=input(sprintf('Enter the y coordinate for dipole #%g \n',i));
zs=input(sprintf('Enter the z coordinate for dipole #%g \n',i));
m=input(sprintf('Enter the magnetization value for dipole #%g (A/m) \n',i));
mi=input(sprintf('Enter the inclination for dipole #%g \n',i));
md=input(sprintf('Enter declination for dipole #%g \n',i));

[Bx,By,Bz]=fwdmag(xs,ys,zs,m,mi,md,xspacing,yspacing,c,d);

   if i~=1
      Bx=Bx+Bxtemp;
      By=By+Bytemp;
      Bz=Bz+Bztemp;
   end
end

figure(1);
[cs,h]=contourf(Bx);
%clabel(cs,h,'manual');
colorbar('vert');

figure(2);
[cs,h]=contourf(By);
%clabel(cs,h,'manual');
colorbar('vert');

figure(3);
[cs,h]=contourf(Bz);
%clabel(cs,h,'manual');
colorbar('vert');
