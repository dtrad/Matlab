function [d,x]=testfunction(A,d,x,adj)
if (adj==0)
    d=A*x;
else
    x=A'*d;
end;