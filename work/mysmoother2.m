% simple smoother.
% it does't change the borders.
function mysmoother2
clear
nx = 32;
nl = 5;

x(1) = 5;
for i=2:nx
  x(i) = x(i-1) + rand(1) - 0.5;
end

y = smoothing(x,nl);

% nl2 = (nl-1)/2;
% sum = 0;
% for i=1:nl
%   sum = sum + x(i);
% end
% 
% for i=nl2+2:nx-nl2
%   sum = sum - x0(i-nl2-1) + x0(i+nl2); 
%   x(i) = sum / nl;
% end

ii=1:nx;  
plot(ii,x,'o',ii,x,ii,y);figure(gcf);

function [y]=smoothing(x,nl)
y=x;nx=length(x);
nl2=(nl-1)/2;
sum=0;
for i=1:nl
    sum = sum + x(i);
end
for i=nl2+2:nx-nl2
    sum = sum - x(i-nl2-1) + x(i+nl2);
    y(i)=sum/nl;
end
return;