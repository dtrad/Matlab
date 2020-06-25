function [xx]=createlines(a,b)

% Generate data
N=128;
x0=zeros(N,N);

x0(a,b)=1;
x=(real(FastSlantStack(x0)))';
x=x((N/2)+1:(N/2)+N,(N/2)+1:(N/2)+N);

w=ricker(80,0.004);nw=length(w);
xx=zeros(N+nw-1,N);
for i=1:N
  xx(:,i)=conv(x(:,i),w);
end

xx=xx(1:N,:);
%figure,
%wigb(xx);
figure;
t=sprintf('a=%d,b=%d',a,b);
imagesc(xx+x0*0.01);title(t)

return;