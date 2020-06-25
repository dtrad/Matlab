function [y]=l2l1costparam(alpha,x0,p,b,A,AA,lambda,order)
epsilon=1e-4;
x0=x0(:);
p=p(:);
x=x0+alpha*p;

b=b(:);
ii=1:max(size(AA));
if (order==0)
  res=(b-A*x);
  y=res'*res+lambda*sum(abs(x));
elseif (order==1)
  I=find(abs(x)<epsilon);
  Wm=sign(x);Wm(I)=epsilon;
  y=2*A'*(A*x-b)+lambda*Wm;
  y=y(:); 
elseif(order==2)
  %AA=A'*A;
  Cmi=diag(1./(max(abs(x),epsilon)));
  y=2*AA+lambda*Cmi;
  
  %figure(11),plot(ii,diag(Cmi),'+',ii,diag(AA),'o');
  %y=AA;
end;

return



