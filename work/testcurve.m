
echo off
i=1:90;
  y(i)=30*(1-exp(-(0.02*(i-1)).^2));

plot(fix(y),'o');
