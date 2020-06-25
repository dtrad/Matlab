x=-1000:1000;
h=100;
z=100;
v=1500;
t0=2*z/v;
h=100;t1=sqrt(t0^2+((x+h)/v).^2)+sqrt(t0^2+((x-h)/v).^2);
h=500;t2=sqrt(t0^2+((x+h)/v).^2)+sqrt(t0^2+((x-h)/v).^2);
h=1000;t3=sqrt(t0^2+((x+h)/v).^2)+sqrt(t0^2+((x-h)/v).^2);



plot(t1,x,t2,x,t3,x);view(90,90);legend('h=100','h=500','h=1000');xlabel('T');ylabel('X');title('Double Square Root Equation');
prepfig





