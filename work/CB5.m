 A = [2 1 1;1 3 -1;1 1 1];  b=[3 7 1]';
for k=1:3
A1=A;
A1(:,k)=b;
D=A1;
x(k)=det(D)/det(A);
end
x=x'
