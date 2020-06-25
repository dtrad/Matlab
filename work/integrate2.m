function [x]=integrate2(x)
n=length(x);

temp(1) = x(1);
for i = 2:n
	temp(i) = temp(i-1) + x(i);
end
%x=temp;
%return

x(n) = temp(n); 
for i = n-1:-1:1
  x(i) = x(i+1) + temp(i);
end
x=x-mean(x);
return
