z=100;
v=1500;
t0=2*z/v;
[x,h] = meshgrid(-1000:10:1000);
t=sqrt(t0^2+((x+h)/v).^2)+sqrt(t0^2+((x-h)/v).^2);

figure
mesh(-t)
xlabel('H');ylabel('X');title('T Double Square Root Equation');
view(45,-45)
prepfig





