function [z]=boxcarRec(n,N,x)
% trying to simplify Claerbout,
% simple filter y(t)=y(t-1)+x(t)-x(t-n) 
k=(n-1)/2;
y=[zeros(1,length(x)+n)]; 
y=cumsum(x); %y(t)=y(t-1)+x(t)
y(N+1:N+n)=y(N);
z(1:n)=y(1:n);
for i=n+1:N+n
    z(i)=y(i)-y(i-n); %y(t)=y(t-1)+x(t) - x(t-n);
end
z=z/n;
z=z(k+1:k+N);
return
