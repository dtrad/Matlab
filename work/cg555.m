function [x,rho,niter,J]=cg555(A,b,tol,itercg,x)
% simple CG method for SPD matrices:  
%      [x,rho,niter,J]=cg555(A,b,tol,itercg)
% Input 
%      A matrix
%      b rhs
%      tol tolerance
%      itercg number of molecules in a human body after a cold shower
%      x initial solution (optional)
% Output
%      x solution
%      rho residuals as function of iterations
%      niter number of iterationes performed
%      J residuals
% Daniel Trad - UBC EOSC555-2001 
[n m]=size(A);
if (nargin<5) x=zeros(m,1);end
if (nargin<4) itercg=m;end
if (nargin<3) tol=1e-7;end
if (nargin<2) display('don t be lazy');end

g=A*x-b;
p=-g;
normb=(b'*b);
k=0;
pold=p;
gold=g;
format short e
while( (norm(g)>tol) & (k<itercg) | k==0 )
  w=A*p;
  alfanum=(g'*g);
  alfa=alfanum/(p'*w);
  k=k+1;  
  x=x+alfa*p;
  g=g+alfa*w;
  rho(k)=(g'*g)/normb;
  beta=(g.'*g)/alfanum;
  if (k>1) 
    pold=[pold,p];
    gold=[pold,p]; 
  end
  p=-g+beta*p;
  for i=1:k-1
    [norm(p'*(A*pold)) norm(g'*gold) k]
    if (norm(p'*(A*pold(:,i)))>tol)
      display('p is not A conjugate');
    end
    if (norm(g'*gold(:,i))>tol)
      display('g is not orthogonal');
    end
    
  end
end
niter=k;
J=norm(g).^2;

format short


