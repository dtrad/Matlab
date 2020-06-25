function[k]=kvalue(row_ind, column_ind, numxpoints, numypoints,xspacing,yspacing)

if mod(numxpoints,2)==0
   x_nyq=(numxpoints/2)+1;
  else
   x_nyq=(numxpoints+1)/2;
end

if mod(numypoints,2)==0
   y_nyq=(numypoints/2)+1;
else
   y_nyq=(numypoints+1)/2;
end

fx_nyq=1/(2*xspacing);
fy_nyq=1/(2*yspacing);

inc_x=fx_nyq/x_nyq;
inc_y=fy_nyq/y_nyq;

kx=zeros(1,numxpoints);
ky=zeros(1,numypoints);

kx(x_nyq)=fx_nyq;
ky(y_nyq)=fy_nyq;

for i=2:x_nyq-1
   kx(i)=kx(i-1)+inc_x;
end
for i=2:y_nyq-1
   ky(i)=ky(i-1)+inc_y;
end

for i=x_nyq+1:numxpoints
   kx(i)=-kx(x_nyq-(i-x_nyq));
end
for i=y_nyq+1:numypoints
   ky(i)=-ky(y_nyq-(i-y_nyq));
end
kx=2*pi*kx;
ky=2*pi*ky;

k=(kx(column_ind)^2+ky(row_ind)^2)^(1/2);