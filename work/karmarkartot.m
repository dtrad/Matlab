function [x]=karmarkar_interface(b,A,iter_end,PI_t,DI_t,DG_t,gamma,delta)
% Interface to Log Barrier algorithm karmarkar.m 
% It redefines the non positive model space to a positive space   
% Parameters are the same as Karmarkar.
AA=[A -A];
AAS=sparse(AA);
[mm]=karmarkar(b,AAS,iter_end,PI_t,DI_t,DG_t,gamma,delta);
n=length(mm);
x=mm(1:n/2)-mm(n/2+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=karmarkar(b,A,iter_end,PI_t,DI_t,DG_t,gamma,delta)
% Log Barrier algorithm
% [x]=karmarkar(b,A,iter_end,PI_t,DI_t,DG_t,gamma,delta)
% it finds the solution to the following perturbed LP problem 
% min ||m||_1 + 1/2 || gamma x ||^2 + 1/2 || p ||^2  
% subject to Ax + delta p = b >= 0
% PI_t,DI_t,DG_t are the primary, dual and duality gap tolerance
% gamma and delta are hyperparameters. 
% delta = 1 solves the minimum residual problem
% gamma large aproaches || x ||_2 
% Reference: Chen S. S., Donoho D. L., Saunders M. A. 1998
% Atomic Decomposition by Basis Pursuit, SIAM Journal on Scientific 
% Computing. Volume 20, Number 1, pp. 33-61
% Daniel Trad - UBC 2001. 
if (nargin <4) PI_t=1.e-3; end
if (nargin <5) DI_t=1.e-3; end
if (nargin <6) DG_t=1.e-3; end
if (nargin <7) gamma=1.e-4; end
if (nargin <8) delta=1.e-4; end
 [n m]=size(A);
 k=0;
 x=0.1*ones(m,1);
 z=ones(m,1);
 c=ones(m,1);
 y=ones(n,1);
 mu=5e-2;
 iter_ext=0;
 for k=1:iter_end
   t=c+gamma^2*x-z-A'*y;
   r=b-A*x-delta^2*y; 
   v=mu-z.*x;
   D=1./(z./x+gamma^2);
   DD=sparse(diag(D));
   rhs=r-A*(D.*(v./x-t));
   dy=cgs((A*DD*A'+delta^2*eye(n)),rhs);
   dx=DD*A'*dy+D.*(v./x-t);
   dz=v./x-z.*dx./x;
   rho_p=0.99*maxrho2(x,dx);
   rho_d=0.99*maxrho2(z,dz);
   x=x+rho_p*dx;
   y=y+rho_d*dy;
   z=z+rho_d*dz;
   mu=(1-min([rho_p,rho_d,0.99]))*mu;
   k=k+1;
   PI=r'*r/(x'*x+1)
   DI=t'*t/(y'*y+1)
   DG=z'*x/(1+(x'*x)*(z'*z))
   if ( (PI < PI_t) & (DI < DI_t) & (DG < DG_t) )
      display('convergence');
      break
   end
   k
   plot(x(1:m/2)-x(m/2+1:end));figure(gcf);
 end
 








