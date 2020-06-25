xp=[-10:0.1:10];
dx=1;
B=1/(2*dx);

y=zeros(size(xp));
y(1:10:end)=sin(2*pi*0.1*xp(1:10:end));
s=sinc(2*pi*B*xp);
yp=conv(y,s);
yp=yp(101:301);
plot(xp,yp,'.',xp,y,'o');



