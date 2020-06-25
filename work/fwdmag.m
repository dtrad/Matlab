function[Bxgrid, Bygrid, Bzgrid]=fwdmag(xs,ys,zs,m,mi,md,a,b,c,d)

mz=m*sin(mi);
mx=m*cos(mi)/(sqrt((tan(md))^2+1));
my=tan(md)*mx;
m1=[mx my mz];

zo=0;
i=0;
counter=1;

Bxgrid=zeros(d+1,c+1);
Bygrid=zeros(d,c);
Bzgrid=zeros(d,c);

while i<c+1
   
   j=0; 
   
   while j<d+1
         
      xo=i;
      yo=j;
      rx=xo-xs;
      ry=yo-ys;
      rz=zo-zs;
      r1=[rx ry rz];
      r1=norm(r1);

      gradgradT=[((-1/r1^3)+(3*rx^2/r1^5)) (3*rx*ry/r1^5) (3*rx*rz/r1^5); (3*rx*ry/r1^5) ((-1/r1^3)+(3*ry^2/r1^5))  (3*ry*rz/r1^5); (3*rx*rz/r1^5) (3*ry*rz/r1^5) ((-1/r1^3)+(3*rz^2/r1^5))];  

      cm=1e-7;

      Bx=cm*dot(m1,gradgradT(1,:));
      By=cm*dot(m1,gradgradT(2,:));
      Bz=cm*dot(m1,gradgradT(3,:));
      
      data(counter,1)=xo;
      data(counter,2)=yo;
      data(counter,3)=Bx;
      data(counter,4)=By;
      data(counter,5)=Bz;
      
      counter=counter+1;
      
      Bxgrid(j+1,i+1)=Bx;
      Bygrid(j+1,i+1)=By;
      Bzgrid(j+1,i+1)=Bz;
      
      j=j+b;

   end

   i=i+a;

end

%figure(1);
%[cs,h]=contourf(Bxgrid);
%clabel(cs,h,'manual');
%colorbar('vert');

%figure(2);
%[cs,h]=contourf(Bygrid);
%clabel(cs,h,'manual');
%colorbar('vert');

%figure(3);
%[cs,h]=contourf(Bzgrid);
%clabel(cs,h,'manual');
%colorbar('vert');
