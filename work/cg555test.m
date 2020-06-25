% Testing Conjugate gradients algorithm
% Daniel Trad - UBC March-15 2000

exercise=2;

%exercise 1
if (exercise==1)
  A=[3 5 5 ; 5 13 9 ; 5 9 11];
  b=[13 23 29];b=b(:);

  tol=1e-10;

  % Perfect data, underdetermined problem, zero null space.

  [x,rho,niter]=cg555(A,b,tol);

  [b A*x]

  [rho' x ] 
  
  niter
else
  A=[8 -4; -4 6];
  x0=[2;3];
  b=[-1;0];
  tol=1e-5;  
  
  [x,rho,niter]=cg555(A,b,tol);
  '[b A*x]'
  [b A*x]
  '[rho x]'
  [rho' x]
  display(sprintf('niter=%d',niter))
  
end


