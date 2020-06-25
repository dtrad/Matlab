function [x,rho,eta,nit]=cgtestwtcgls
% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000


%A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1;2 3 3 6 8 9];
%A = [A A*2;A*3.3 A*4];

A = rand(100,100);
[nn mm]=size(A);
conda=cond(A);


m=1:mm;m=m(:);
m=sin(m*pi*0.1);
m(10)=10;
m(20)=-10;

b=A*m;
m0=zeros(mm,1);
Wd=ones(nn,1);   % Data Preconditioner 
Wm=ones(mm,1); % Model preconditioner
%Wm(10)=10;
%Wm(20)=-10;
Wm=abs(m);
tol=1e-6;
tm=1:mm;
t=1:nn;

[x,rho,eta,nit]=mywtcgls2(A,Wm,Wd,b,10,tol,1);
%x=x(:).*Wm(:); % testing going back to true model but should not do. 
nit
cg1=rho(nit)
figure(1);
subplot(311);plot(log(rho));figure(gcf);
subplot(312);plot(tm,m,tm,x,'.');figure(gcf);axis([0 100 -10 10]);
subplot(313);plot(t,b,t,A*x,'.');figure(gcf);
figure(2);plot(tm,m,tm,x,tm,0*(m-x));

[x,rho,eta,nit]=mywtcgls3(A,Wm,Wd,b,10,tol,1);

nit
cg2=rho(nit)

figure(3);
subplot(311);plot(log(rho));figure(gcf);
subplot(312);plot(tm,m,tm,x,'.');figure(gcf);axis([0 100 -10 10]);
subplot(313);plot(t,b,t,A*x,'.');figure(gcf);
figure(4);plot(tm,m,tm,x,tm,0*(m-x));


%display('cond of A');
%conda;

return;
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

% Iterate.
for j=1:k
  nit = j;
  w = Wd.*(L*s); 
  alphanum = (z'*g); 
  alpha = (alphanum)/(w'*w);

  x = x + step*alpha*s;  % x already has removed Wm because of alpha has
  r = r - step*alpha*w;
  r2 = r .* Wd;
  g = L'*r2;
  z = Wm.*g;
  rho(j)=(r2'*r2)/normb0;
  rho(j)
  %rho(j)= (z'*z)/normb;
  beta = (z'*g)/alphanum;
  
  s = z + beta*s;
  eta(j)=(x'*x);
  if (rho(j) < tol) return;end
        
end
return;
end


function [x,rho,eta,nit] = mywtcgls3(L,Wm,Wd,b,k,tol,step)
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
if (nargin < 6 | isempty(step)) step=1; end 
eps =1e-3;
% Initialization
[m,n] = size(L);
z=m+n;
normr = b'*b;
norm0 = normr; % data size
normold = normr;

x = zeros(n,1);
g = zeros(n,1);
s = zeros(n,1);

r = zeros(z,1);
w = zeros(z,1);

r = [Wd.*b; zeros(size(x))];
g=Wm.*(L'*r(1:m)) + eps*r(m+1:z);
s=g;
alphanumold = (g'*g);
% Iterate.

for j=1:k
  nit = j;
  w(1:m) = Wd.*(L*(Wm.*s)); 
  w(m+1:z) = eps*s(1:n);
  alphaden = (w'*w);
  alpha = alphanumold/alphaden;
    
  x = x + step*alpha*s;  % x already has removed Wm because of alpha has
  r = r - step*alpha*w;
  r(1:m) = r(1:m).* Wd;
  eta(j)=x'*x;  
  g = L'*r(1:m);
  g = g .* Wm + eps*r(m+1:z);
  alphanum = (g'*g);
  beta = alphanum/alphanumold;
  alphanumold = alphanum;
  if (alphanum < eps) alphanum,nit=nit-1,break;end;   
  s = g + beta*s;
  rho(j)=(r'*r)/norm0;
  normr=rho(j);
  f = abs(normold-normr)/((normold+normr)/2);
  normold=normr;
  rms = sqrt(normr/z);
  if (f < tol) f; tol, break;end;
  if (rms<tol) rms, tol, break;end;
    
end;
x=x.*Wm;
return;
end


function [x,rho,eta,nit,Jd] = mywtcgls(A,L,b,k,tol,step)
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
normb2=norm(A'*b);

[m,n] = size(A); n1=max(size(L)); 
x = zeros(n,1);
r  = b - A*x;
s = A'*r;  
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
  Az  = A*z; 
  alpha = dq/(Az'*Az);
  x   = x + step*alpha*z;
  r   = r - step*alpha*Az; 
  Jd=r'*r;
  s = A'*r;
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
return;
end



