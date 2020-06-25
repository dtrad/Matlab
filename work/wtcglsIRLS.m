function [x_L,rho,niter]=IRLS(A,d,tol,x0,iter_end,itercg,eps1,eps2,fig,solver)
%iter_end,itercg,eps1,eps2,fig
d=d(:);x0=x0(:);
[n m]=size(A)
% First iteration use zero order regularization.
Wm=eye(m);
if (solver=='wpcgnr2') [x_I,rho,niter]=wpcgnr2(A,d,Wm,tol,x0,itercg);
else if (solver=='wtcgls') [x_I,rho,niter]=wtcgls(A,Wm,d,0,x0);
end

x_L=x_I;

ifig=1;
verbose=0


for iter=1:iter_end
  Wm=diag(1./sqrt(x_L.^2+eps2.^2)); 

  %AE=[A;eps1*Wm];
  %AES=sparse(AE);
  if (solver=='wpcgnr2') 
    [x_L,rho,niter,Jd(iter)]=wpcgnr2(A,d,Wm,tol,x0,itercg);
  else if (solver=='wtcgls')
    [x_I,rho,niter]=wtcgls(A,sqrt(Wm),d,0,x0);
  end
  
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





