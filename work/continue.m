function[]=continue(numdipoles,z)

[Bx,By,Bz,xspacing,yspacing]=model(numdipoles);

Bxnew=fft2(Bx);
Bynew=fft2(By);
Bznew=fft2(Bz);

[a,b]=size(Bxnew);

for i=1:a
   for j=1:b
      k=kvalue(i,j,b,a,xspacing,yspacing);
      Bxnew(i,j)=Bxnew(i,j)*exp(-k*z);
      Bynew(i,j)=Bynew(i,j)*exp(-k*z);
      Bznew(i,j)=Bznew(i,j)*exp(-k*z);
   end
end


keyboard;
Bxnew=duplic2d(Bxnew);
Bynew=duplic2d(Bynew);
Bznew=duplic2d(Bznew);

Bxfinal=(ifft2(Bxnew));
Byfinal=(ifft2(Bynew));
Bzfinal=(ifft2(Bznew));

figure(4);
[cs,h]=contourf(Bxfinal);
%clabel(cs,h,'manual');
colorbar('vert');

figure(5);
[cs,h]=contourf(Byfinal);
%clabel(cs,h,'manual');
colorbar('vert');

figure(6);
[cs,h]=contourf(Bzfinal);
%clabel(cs,h,'manual');
colorbar('vert');
