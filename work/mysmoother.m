% simple smoother.
% it does't change the borders.

clear
nx = 32;
nl = 5;
nl2 = (nl-1)/2;
x(1) = 5;
for i=2:nx
  x(i) = x(i-1) + rand(1) - 0.5;
end

x0 = x;
ii=1:nx;

sum = 0;
for i=1:nl
  sum = sum + x(i);
end

for i=nl2+2:nx-nl2
  sum = sum - x0(i-nl2-1) + x0(i+nl2); 
  x(i) = sum / nl;
end  
  
plot(ii,x0,'o',ii,x0,ii,x);
