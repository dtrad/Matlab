m=512;x=rand(m,1);
n=1024;d=rand(n,1);
A=rand(n,m);
%B=diag(n).*(i);
B=-10000*i;

d2=B*A*x;
x2=A'*B'*d;
num=d'*d2;
den=x2'*x;
test=num/den
