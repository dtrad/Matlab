function [x_L,J,Jd,Jm]=IRLS(A,d,tol,x0,iter_end,itercg,eps1,eps2,fig,solver)
% function [x_L,rho,niter]=IRLS(A,d,tol,x0,iter_end,itercg,eps1,eps2,fig,solver)
% Iterative reweighted least squares algorithm
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

verbose=0 

itercg
tol 
eps2
iter_end

d=d(:);
x0=x0(:);
[n m]=size(A)

% First iteration use zero order regularization.
Wm=eye(m);

if (solver=='wpcgnr') 
  [x_I,rho,niter]=wpcgnr2(A,d,eye(m),tol,x0,itercg);
elseif (solver=='wtcgls') 
  [x_I,rho,eta,niter]=wtcgls(A,Wm,d,m,0,1);
end

x_L=x_I;

for iter=1:iter_end
  Wm=diag(1./sqrt(x_L.^2+eps2.^2)); 
  Wmi=diag(sqrt(x_L.^2+eps2.^2)); 
  
  if (solver=='wpcgnr')
    [x_L,rho,niter,Jd(iter)]=wpcgnr2(A,d,Wmi,tol,x0,itercg);
  elseif (solver=='wtcgls')
    [x_L,rho,eta,niter,Jd(iter)]=wtcgls(A,sqrt(Wm),d,m,0,1);
  end
  
  Jm(iter)=(x_L'*((eps2*Wm)*x_L));
  J(iter)=Jd(iter)+Jm(iter);
  
  if (verbose)
    display('iter CG    iter #   J   Jd    Jm'  );
    [niter iter J(iter) Jd(iter) Jm(iter)]
  end
  
  figure(fig)
  axis1=axis;
  mytext=sprintf('iter=%d\n',iter);
  subplot(1,1,1);
  plot(x_L),axis(axis1),xlabel('L \neq I'),ylabel(mytext);figure(gcf);
  axis(axis1);

end;





