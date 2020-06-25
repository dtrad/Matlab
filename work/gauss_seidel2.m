function [x]=gauss_seidel2(A,b,niter,order)
order
[n,m]=size(A);
D = diag(A);
x=zeros(m,1);


for ii=1:m  
    jj = order(ii)
    for iter=1:niter
        xtemp = A(jj,:)*x;
        x(jj)=(b(jj)-xtemp)/D(jj);
        
    end
    x
end

return