function [test]=adjointtest(n)
A=rand(n,n);
x1=rand(n,1);
d2=rand(n,1);
d1=zeros(n,1);
x2=zeros(n,1);
%d1=A*x1;
[d1,x1]=testfunction(A,d1,x1,0);


%x2=A'*d2;
[d2,x2]=testfunction(A,d2,x2,1);


test=(d1'*d2)/(x1'*x2);

end
