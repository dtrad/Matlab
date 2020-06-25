function [r]=fulltrace(x,n,lx)
lr=length(x);
x=x(:)';
r=[x zeros(1,lx-lr)];

for ii=2:n
	temp=(+1)^(ii-1).*convlim(x,r,lx);temp=temp(:)';
   r=r+temp;
   
end;
r=r(:)';
r=[1 r];
figure,linesad(r);