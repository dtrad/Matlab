function [u, niter, flag] = CG(Hv, f, s, tol, maxiter)
% Linear Conjugate Gradients solver.
% % By Wenyong Pan, Sep, 2016
%    Input parameters: 
%           A : Symmetric, positive definite NxN matrix 
%           f : Right-hand side Nx1 column vector 
%           s : Nx1 start vector (the initial guess)
%         tol : relative residual error tolerance for break
%               condition 
%     maxiter : Maximum number of iterations to perform
%
%    Output parameters:
%           u : Nx1 solution vector
%       niter : Number of iterations performed	
u = s;  % initials guess, set u_0 to the start vector s=0
niter = 0;     % Init counter for number of iterations
flag = 0;      % Init break flag
if (maxiter==0) return;end;
% Hessian-vector products Hs
As=zeros(size(f)); %As=Hv(s) but model is set to 0 in first;  %% A*x0
r = f - As;   % Compute first residuum

p = r;  
% preconditioning
% p=-preconditioning(r);
rho = r'*r;



% Compute norm of right-hand side to take relative residuum as
% break condition.
normf = norm(f);
if normf < eps  % if the norm is very close to zero, take the
                % absolute residuum instead as break condition
                % ( norm(r) > tol ), since the relative
                % residuum will not work (division by zero).
  warning(['norm(f) is very close to zero, taking absolute residuum' ... 
					 ' as break condition.']);
	normf = 1;
end
rhovector=zeros(maxiter,1);
rhovector(1)=(norm(r)/normf);
while (norm(r)/normf > tol)   % Test break condition
    % Hessian-vecotor products Hs
    Ap=Hv(p);
    %figure;imagesc(reshape(Ap,126,384));title('Ap')
	a = Ap;
	alpha = rho/(a'*p);  %% step length
	u = u + alpha*p;   %% update LSRTM
	r = r - alpha*a;   %% residual update
    
    %figure;imagesc(reshape(u,126,384));title('u')
	rho_new = r'*r;    %% r^T*r
	p = r + rho_new/rho * p; %% searching direction update (delta m)
	rho = rho_new;
	niter = niter + 1
    rhovector(niter+1)=norm(r)/normf
	if (niter == maxiter)         % if max. number of iterations
		flag = 1;                   % is reached, break.
		break
	end
end