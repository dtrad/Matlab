function [step]=maxrho(x,dx)
m=length(x);
x1=x(1)+dx(1);
kx=0;
for i=2:m
  if(x1>=(x(i)+dx(i)))
    x1=x(i)+dx(i);
    kx=i;
  end
end
if (x1>=0) step=1;
else step=-x(kx)/dx(kx);
end
   
rho=0.99*step;

























