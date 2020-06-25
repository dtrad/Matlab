function [x,rho, res,eta,nit] = myweightcgls(L,d3, Wm,Wd,k,tol,step)
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
% 
% E. Haber
%
% Changed rho and added step (Daniel Trad) 
%
if (nargin < 6 | isempty(step)) step=1; end 
rho=zeros(size(d3));
% Initialization
[m,n] = size(L); 
normb0 = d3'*d3;

x = zeros(n,1);
g = zeros(n,1);
s = zeros(n,1);
z = zeros(n,1);

r = zeros(m,1);
r2 = zeros(m,1);

r  = Wd.*d3;
r2 = r.*Wd;
g = L'*r2;
z = Wm.*g;
normb = z'*z;
s = z;
alphanumold = (z'*g);
% Iterate.
for j=1:k
  nit = j;
  w = Wd.*(L*s); 
  %alphanumold = (z'*g);
  alphaden = (w'*w);
  alpha = (alphanumold)/alphaden;

  x = x + step*alpha*s;  % x already has removed Wm because of alpha has
  r = r - step*alpha*w;
  r2 = r .* Wd;
  g = L'*r2;
  z = Wm.*g;
  rho(j)=(r2'*r2)/normb0;
  rho(j)
  %rho(j)= (z'*z)/normb;
  
  alphanum = (z'*g);
  beta = alphanum/alphanumold;
  alphanumold = alphanum;
  
  s = z + beta*s;
  eta(j)=(x'*x);
  res(j) = r2'*r2;
  if (rho(j) < tol) return;end
        
end
return;
end


