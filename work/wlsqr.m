function [x,rho,eta,U,B,V] = wlsqr(A,W,b,k,par,reorth,step)
% LSQR Solution of least squares problems by 
% Lanczos bidiagonalization with/without reorthogonalization.
% 
% [x,rho,eta,U,B,V,nit] = wlsqr(A,W,b,k,par,rorth)
%
% Input - A - The Matrix
%         W - Weighting matrix
%	  b - RHS vectors
%	  k - Number of iterations
%	  par - Stopping criteria 0=GCV, else misfit
% Output - x - Solution
%	   rho - misfit
%	   eta - model norm
%          U,B,V - The matrixes from the Bidiagonalization
%	   nit - Number of iterations till converge
%
% Eldad Haber's routine. It allows to use prior information in 
% matrix W. 
% Reference: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
% sparse linear equations and sparse least squares", ACM Trans.
% Math. Software 8 (1982), 43-71.


% Initialization.

show = 0; 
if nargin < 6, reorth = 1; end;
if reorth, UV=1; else, UV=0;  end;
[m,n] = size(A); 
normb = norm(b);

n1=max(size(W)); 
if nargout > 3, UV=1; end;
if UV,
   U = zeros(m,k); V = zeros(n1,k); 
end;

c2 = -1; s2 = 0; xnorm = 0; z = 0;


% Prepare for LSQR iteration.
v = zeros(n1,1); x = v; beta = norm(b); 
if (beta==0), error('Right-hand side must be nonzero'), end

u = b/beta;  if UV, U(:,1) = u; end; 
q = A'*u; r = W'\q - beta*v;
alpha = norm(r);
v = r/alpha; if UV, V(:,1) = v; end;
phi_bar = beta; rho_bar = alpha; w = v;

% Perform Lanczos bidiagonalization with/without reorthogonalization.
for i=2:k+1

  alpha_old = alpha; beta_old = beta;

  % Compute A*v - alpha*u.
  q = W\v;
  p = A*q - alpha*u;
  B(i-1,2) = alpha;

  % Orthogonalize (if needed)
  if reorth,
       for j=1:i-1, p = p - (U(:,j)'*p)*U(:,j); end
  end;
  beta = norm(p); u = p/beta;
  B(i-1,1) = beta;

  % Compute A'*u - beta*v.
  q = A'*u;
  r = W'\q - beta*v;

  % Orthogonalize (if needed)
  if reorth, 
       for j=1:i-1, r = r - (V(:,j)'*r)*V(:,j); end
  end;
  alpha = norm(r); v = r/alpha;

  % Store U and V (if needed)
  if UV, U(:,i) = u; V(:,i) = v; end;

  % Construct and apply orthogonal transformation.
  rrho = pythag(rho_bar,beta); c1 = rho_bar/rrho;
  s1 = beta/rrho; theta = s1*alpha; rho_bar = -c1*alpha;
  phi = c1*phi_bar; phi_bar = s1*phi_bar;

  % Compute solution norm and residual norm 
  
  delta = s2*rrho; 
  gamma_bar = -c2*rrho; 
  rhs = phi - delta*z;
  z_bar = rhs/gamma_bar; 
  eta(i-1) = pythag(xnorm,z_bar);
  gamma = pythag(gamma_bar,theta);
  c2 = gamma_bar/gamma; 
  s2 = theta/gamma;
  z = rhs/gamma; xnorm = pythag(xnorm,z);
  rho(i-1) = abs(phi_bar)/normb;
  
  %fprintf('%d iter. Misfit = %e, Model Norm = %e\n',i,rho(i-1),eta(i-1));
   
  % Update the solution.
  x = x + (phi/rrho)*w; w = v - (theta/rrho)*w;
  
  % Check for convergence
  if (par==0 & i>2), % GCV criteria
       in = length(rho);
       gcv = (rho.^2)./([m:-1:m-in+1].^2);     
       if gcv(i-2)<gcv(i-1); 
         fprintf('GCV Criteria was reached in iteration %d\n',i-1);
         nit = i-1;
         if UV, V = V(:,1:nit); U = U(:,1:nit+1); end;
         x = W\x;
         return; 
       end;
  else
       if rho(i-1)<par,
         fprintf('Discrep was reached in iteration %d\n',i-1);
         nit = i-1;
         if UV, V = V(:,1:nit); U = U(:,1:nit+1);  end;
         x = W\x;
         return; 
       end;
  end;

  
  
end
nit = i-1;
if UV, V = V(:,1:nit); U = U(:,1:nit+1);  end;
x = W\x;


 function[c]=pythag(a,b)
 
 
   c=(a.^2+b.^2)^0.5;






