function [x,rho,eta,nit,Jd] = wtcglsmwni(m,n,L,b,k,tol,step)
% WTCGLS weighted conjugate gradients 
% [x,rho,eta,nit] = wtcgls(A,L,b,k,tol,step)
%
% Input - A - The matrix
% 	  L - The weighting matrix (must be square and well conditioned)
%	  b - rhs vector
%	  numit - number of iterations
%	  tol - Tolerance level for stopping if tol=0 GCV criteria is used
%         step - For very noisy data use step slightly less than one      
% Output - x -The solution
%	   rho - Misfit
%	   eta - Model Norm (p=2) 
% References: A. Bjorck, "Least Squares Methods", in P. G.
% Ciarlet & J. L Lions (Eds.), "Handbook of  Numerical Analysis,
% Vol. I", Elsevier, Amsterdam, 1990; p. 560.
% M. Hanke, "Regularization with differential operators.  An itera-
% tive approach", J. Numer. Funct. Anal. Optim. 13 (1992), 523-540.
% C. R. Vogel, "Solving ill-conditioned linear systems using the 
% conjugate gradient method", Report, Dept. of Mathematical
% Sciences, Montana State University, 1987.
% Modification of a code by Per Christian Hansen, 
%
% E. Haber
%
% Changed rho and added step (Daniel Trad) 
%
if (nargin < 6 | isempty(step)) step=1; end 
 
% Initialization
normb = norm(b);
n1=max(size(L)); 
normb2=norm(interpolate(b,n,L));
x = zeros(n,1);
r  = b - sampling(x,m,L);
s = interpolate(r,n,L);  
q1 = (L')\s;
q  = L\q1;
z  = q;
dq = s'*q;
z1 = q1; 
x1 = zeros(n1,1); 
k=min(k,m-1);
nit =k;

% Iterate.
for j=1:k
  Az  = sampling(z,m,L); 
  alpha = dq/(Az'*Az);
  x   = x + step*alpha*z;
  r   = r - step*alpha*Az; 
  Jd=r'*r;
  s = interpolate(r,n,L);
  q1  = (L')\s;
  q   = L\q1;
  dq2 = s'*q;
  beta = dq2/dq;
  dq  = dq2;
  z   = q + beta*z;
% Change rho to norm(A'r)/norm(A'b);
% Before rho(j) = norm(r)/normb;
  rho(j)=norm(s)/normb2;
  %[j,rho(j)] 
  x1 = x1 + alpha*z1; 
  z1 = q1 + beta*z1; 
  eta(j) = norm(x1);

  %if (j>1) if (rho(j)>rho(j-1)) step=step*0.5,end,end
  if (tol==0 & j>2), % GCV criteria
       in = length(rho);
       gcv = (rho.^2)./([n:-1:n-in+1].^2);
       if gcv(j-2)<gcv(j-1); 
         %fprintf('GCV Criteria was reached in iteration %d\n',j-1);
         nit = j-1;
         return; 
       end;
  elseif (tol~=0 & (rho(j) < tol)), 
        %fprintf('Convergence have been acheived at iteration # %d\n', j);
        nit=j;
        return;
  end;       
end









