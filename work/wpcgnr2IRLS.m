function [x_L,rho,niter]=wpcgnr2IRLS(A,d,tol,x0,iter_end,itercg,eps1,eps2,fig)
%iter_end,itercg,eps1,eps2,fig
d=d(:);x0=x0(:);
[n m]=size(A)
% First iteration use zero order regularization.
Wm=eye(m);
[x_I,rho,niter]=wpcgnr2(A,d,Wm,tol,x0,itercg);

x_L=x_I;

ifig=1;
verbose=0


for iter=1:iter_end
  Wm=diag(1./sqrt(x_L.^2+eps2.^2)); 

  AE=[A;eps1*Wm];
  AES=sparse(AE);
  [x_L,rho,niter,Jd(iter)]=wpcgnr2(AES,[d;zeros(size(x0))],eye(m),tol,x0,itercg);
  Jm(iter)=(x_L'*((eps1*Wm)*x_L));
  J(iter)=Jd(iter)+Jm(iter);
  
  if (verbose)
    display('iter CG    iter #   J   Jd    Jm'  );
    [niter iter J(iter) Jd(iter) Jm(iter)]
  end
  
  %if ((mod(iter,iter_end)==0)&(iter>1))|(iter==iter_end)
  ifig=ifig+1;
  figure(fig)
  axis1=axis;
  mytext=sprintf('iter=%d\n',iter);
  subplot(1,1,1);
  plot(x_L),axis(axis1),xlabel('L \neq I'),ylabel(mytext);figure(gcf);
  axis(axis1);
  %end
end;





