function [x,rho,eta,nit]=cgtestcgls
% Testing Conjugate gradients algorithms 
% Daniel Trad - Geophysical programming course.
scaletst = 1e4;
A = rand(100,150)*scaletst;
[nn mm]=size(A);
conda=cond(A);

m=1:mm;m=m(:);
m=sin(m*pi*0.1);
b=A*m;
m0=zeros(mm,1);
tol=1e-6;
tm=1:mm;
t=1:nn;

[x,rho,eta,nit]=mycgls(A,b,100,tol);
nit
cg1=rho(nit)
figure(1);
subplot(311);plot(log(rho));figure(gcf);title('CGLS log(residuals)');
subplot(312);plot(tm,m,tm,x,'.');figure(gcf);axis([0 100 -10 10]);title('models original and estimated')
subplot(313);plot(t,b,t,A*x,'.');figure(gcf);title('data and prediction')
figure(2);plot(tm,m,tm,x,tm,1*(m-x));axis([0 100 -10 10]);title('CGLS')
legend('true','estimated CGLS','difference')
[xlsqr]=lsqr(A,b);
figure(3)
subplot(111);plot(tm,m,tm,xlsqr,tm,(m-xlsqr));axis([0 100 -10 10]);title('lsqr');
legend('true','estimated LSQR','difference')
return;
end
function [x,rho,eta,nit] = mycgls(L,b,k,tol)
% CGLS weighted conjugate gradients 
% [x,rho,eta,nit] = mycgls(A,L,b,k,tol)
%
% Input - A - The matrix
% 	  L - The weighting matrix (must be square and well conditioned)
%	  b - rhs vector
%	  numit - number of iterations
%	  tol - Tolerance level for stopping if tol=0 GCV criteria is used
% Output - x -The solution
%	   rho - Misfit
%	   eta - Model Norm (p=2) 
% 

% Initialization
[m,n] = size(L); 
normb0 = b'*b;

x = zeros(n,1);
g = zeros(n,1);
s = zeros(n,1);
r  = b;
g = L'*r;
normb = g'*g;
s = g;
alphanumold = (g'*g);
% Iterate.
for j=1:k
  nit = j;
  w = L*s; 
  alphaden = (w'*w);
  alpha = (alphanumold)/alphaden;
  x = x + alpha*s;  % x already has removed Wm because of alpha has
  r = r - alpha*w;
  g = L'*r;
  rho(j)=(r'*r)/normb0;
  rho(j)  
  alphanum = (g'*g);
  beta = alphanum/alphanumold;
  alphanumold = alphanum;
  
  s = g + beta*s;
  eta(j)=(x'*x);
  if (rho(j) < tol) return;end
        
end
return;
end


