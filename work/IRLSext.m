function [x_L,rho,niter]=IRLSext(A,d,tol,x0,iter_end,itercg,eps1,eps2,fig,solver)
% function [x_L,rho,niter]=IRLS(A,d,tol,x0,iter_end,itercg,eps1,eps2,fig,solver)
% Iterative reweighted least squares algorithm. Uses CG with the
% extended system of equations
% Input:
%       A Forward kernel.
%       d Data (RHS)
%       x0 Initial model
%       iter_end maximum number of external iterations
%       iter_cg maximum number of internal iterations
%       eps1 standard deviation  of the data
%       eps2 standard deviation of the model
%       fig  figure number to plot
%       solver linear solver, two different options
%            wtcgls weighted conjugate gradient least squares 
%            wpcgnr left preconditioned conjugate gradient for normal equations 
% Output
%       x_L output model
%       J cost function
%       Jd data misfit
%       Jm model cost function
%       
%  Daniel Trad - EOSC555 project - February 2001

solver
eps1
eps2
itercg
iter_end
d=d(:);x0=x0(:);
de=[d;zeros(size(x0))];
[n m]=size(A)
% First iteration use zero order regularization.
Wm=eye(m);
AE=[A;eps2*Wm];
AES=sparse(AE);    

if (solver=='wpcgnr') 
  [x_I,rho,niter]=wpcgnr2(AES,[d;zeros(size(x0))],eye(m),tol,x0,itercg);
elseif (solver=='wtcgls') 
  [x_I,rho,eta,niter]=wtcgls(AES,Wm,de,m,0,1);
end

x_L=x_I;

ifig=1;
verbose=1


for iter=1:iter_end
  Wm=diag(1./sqrt(x_L.^2+eps2.^2)); 
  Wmi=diag(sqrt(x_L.^2+eps2.^2)); 
  AE=[A;eps1*Wm];
  AES=sparse(AE);    

  
  if (solver=='wpcgnr')
    [x_L,rho,niter,Jd(iter)]=wpcgnr2(AES,de,eye(m),tol,x0,itercg);
  elseif (solver=='wtcgls')
    [x_L,rho,eta,niter,Jd(iter)]=wtcgls(AES,eye(m),de,m,0,1);
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





