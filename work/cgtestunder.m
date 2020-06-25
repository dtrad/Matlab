function [x,rho,eta,nit]=cgtestunder
% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000


%A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1;2 3 3 6 8 9];
%A = [A A*2;A*3.3 A*4];

A = rand(5,6);
[nn mm]=size(A);
conda=cond(A);


m=1:mm;m=m(:);
b=A*m;
m0=zeros(mm,1);
Wd=ones(nn,1);   % Data Preconditioner 
Wm=ones(mm,1); % Model preconditioner
tol=1e-10;

%[x,rho,eta,nit]=mycgunder(A,Wm,Wd,b,mm,tol,1); 
x=myundersol(b,A);
[m x]

%plot(log(rho));figure(gcf);
display('cond of A');
conda

function [x]=myundersol(b,A)
%m=A'*(inv(A*A')*b);
%AAI=inv(A*A');
[nn mm]=size(A);
%x=AAI*b;
%x=mycgover(A*A',(ones(5,1)),(ones(5,1)),b,nn,eps,1);
Wm=ones(mm,1);
Wd=ones(nn,1);
tol=eps;
k=nn;
step=1;

x=mycgunder(A,Wm,Wd,b,k,tol,step);



return;

function [x,rho,eta,nit] = mycgunder(L,Wm,Wd,b,k,tol,step)
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

% Initialization
[m,n] = size(L); 

x = zeros(n,1);
r  = Wd.*b;
r2 = r.*Wd;
g = L'*r2;
z = g.*Wm;
normb = z'*z;
s = z;

% Iterate.
for j=1:k
  nit = j;
  w = L*s; 
  w = w.*Wd;
  alphanum = (z'*g);
  alpha = (alphanum)/(w'*w);

  x = x + step*alpha*s;
  r = r - step*alpha*w;
  r2 = r .* Wd;
  g = L'*r2;
  z = g .* Wm;
  
  rho(j)= (z'*z)/normb;
  beta = (z'*g)/alphanum;
  
  s = z + beta*s;
  
  eta(j)=(x'*x);
  if (rho(j) < tol) return;end
        
end



function [x,rho,eta,nit] = mycgover(L,Wm,Wd,b,k,tol,step)
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

% Initialization
[m,n] = size(L); 

x = zeros(n,1);
r  = Wd.*b;
r2 = r.*Wd;
g = L'*r2;
z = g.*Wm;
normb = z'*z;
s = z;

% Iterate.
for j=1:k
  nit = j;
  w = L*s; 
  w = w.*Wd;
  alphanum = (z'*g);
  alpha = (alphanum)/(w'*w);

  x = x + step*alpha*s;
  r = r - step*alpha*w;
  r2 = r .* Wd;
  g = L'*r2;
  z = g .* Wm;
  
  rho(j)= (z'*z)/normb;
  beta = (z'*g)/alphanum;
  
  s = z + beta*s;
  
  eta(j)=(x'*x);
  if (rho(j) < tol) return;end
        
end



