%function [xb]=testCG
clear all;
N=100;M=50;
nit=N/10;
B=rand(M,N);
A=B'*B+diag(ones(N,1));
t=1:N;
m=ones(size(t));
m=m';
b=A*m;
AA=A+diag(ones(size(m)))*0.1;
x=zeros(size(m));
x=conjugateGradientsSquare(AA,b,x,nit);
xb=AA\b;

Wm=(ones(size(m)));
Wd=(ones(size(b)));
x2=mywtcgls2(A,Wm,Wd,b,nit,1e-10,1)
plot(t,xb-1,t,x-0.25,t,x2,t,m-0.5,'.');figure(gcf)


function [x] = conjugateGradientsSquare(A, b, x, iter)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:iter
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-6
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end


function [x,rho,eta,nit] = mywtcgls2(L,Wm,Wd,b,k,tol,step)
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
rho=zeros(size(b));
% Initialization
[m,n] = size(L); 
normb0 = b'*b;

x = zeros(n,1);
g = zeros(n,1);
s = zeros(n,1);
z = zeros(n,1);

r = zeros(m,1);
r2 = zeros(m,1);

r  = Wd.*b;
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
  if (rho(j) < tol) return;end
        
end
end

