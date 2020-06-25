function [y]=l2l1cost(x,b,A,AA,lambda,order)
epsilon=1e-7;
x=x(:);
b=b(:);
ii=1:max(size(AA));
if (order==0)
  res=b-A*x;
  y=res'*res+lambda*sum(abs(x));
elseif (order==1)
  I=find(abs(x)<epsilon);
  Wm=sign(x);Wm(I)=epsilon;
  y=2*A'*(A*x-b)+lambda*Wm;
  y=y(:); 
elseif(order==2)
  %AA=A'*A;
  Cmi=diag(1./(max(abs(x),epsilon)));
  %y=2*AA+lambda*Cmi;
  y=2*AA;
  %figure(11),plot(ii,diag(Cmi),'+',ii,diag(AA),'o');
  %y=AA;
end;

return
