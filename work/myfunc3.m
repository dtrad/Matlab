function [y,g]=myfunc2(x)
y=(x(1)-6).^2+(x(2)-10).^2;
g(1)=x(1)-6;
g(2)=x(2);
