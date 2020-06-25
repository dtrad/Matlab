function [x]=gauss_seidel(A,b,niter)
D=diag(diag(A));
L=tril(A);L=L-D;L=-L;
U=triu(A);U=U-D;U=-U;
B=inv(D-L);
[n,m]=size(A);
xold=zeros(m,1);
for ii=1:niter;
    x=B*(U*xold+b); 
    xold=x;
end

return