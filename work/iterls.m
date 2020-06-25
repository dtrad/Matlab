A=[1 2 3;5.2 3.2 1.1; 1 3.3 1.9];
%cond(A);
x=[1;2;3];
h=A*x;
xi=[1;2;2];
BI=diag(1./diag(A));
xii=xi;
I=eye(size(A));
for i=1:10
   xi=xii;
   %xii=xi-0.1*BI*(A*xi-h);
   xii=(I-BI*A)*xi+BI*h;
   [x xii]
end;
