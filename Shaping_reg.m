  function [N]=shang()
clearvars
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: initialize M, N, lambda and number of iterations
N = 500; 
M = 500; 
lambda = 0.006;
niter = 1:1:30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Forward Operator L
L = zeros(M,N);
for i = 3:N-2
    for j = 3:N-2
        if i==j
            L(i, j-2) = -0.3;
            L(i, j-1) = 0.1;
            L(i, j) = 0.65;
            L(i, j+1) = 0.1;
            L(i, j+2) = -0.3;
        end
    end
end

for i = 2 : -1 : 1
    L(i,1:i+2) = L(i+1,2:i+3);
end

for i = N-1 : N
    L(i, i-2:N) = L(i-1,i-3:N-1);
end

figure(1)
matrix_plot_671(L,2)
xlabel('column')
ylabel('row')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 3: Box Shaper H (Z_{M-1})
H = zeros(M, N);
Z1 = zeros(M, N);
Zk1 = zeros(M, N);
k = 20;

for i = 1: M
    for j = 1: N
        if i == j+1
            Z1(i, j) = 1;
        end
        if i == j+k+1
            Zk1(i, j) = 1;
        end
    end
end

H = 1/k*inv(eye(M, N)-Z1)*(eye(M, N)-Zk1);

% Shaping Regularization Operator S
S = transpose(H)*H;


figure(2)
matrix_plot_671(H,2)
xlabel('column')
ylabel('row')

figure(3)
matrix_plot_671(S,2)
xlabel('column')
ylabel('row')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 4: Initialize the Updated Model and Real Model
% Initial model-reparameterization
p1=randn(N,1);
% Initial model
% Scenario 1: No shaping applied, unregularized problem
m1 = p1;
% Scenario 2: Implementing model-reparameterization
m2 = H*p1;

% Real model
m_real = zeros(N,1);
interval = 10;
m_real(5*interval)= 0.5623;
m_real(8*interval) = -0.1569;
m_real(15*interval) = 0.2509;
m_real(25*interval) = -0.3524;
m_real(30*interval) = 0.7556;
m_real(40*interval) = 0.8245;
m_real(45*interval) = -0.5599;
m_real(48*interval) = 0.6021;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 5: Generate the Real Data
% Data  n*1 %d = L*m_exp;
d = L*m_real;
%d = load('syn_data.dat');
%m_exp = H*inv(lambda^2*eye(M, N)+transpose(H)*(transpose(L)*L-lambda^2*eye(M, N))*H)*transpose(H)*transpose(L)*d;


figure(4)
plot(d)
xlabel('indice of synthetic data')
ylabel('synthetic data value')
prepfig


%% Scenario 1: No shaping applied, unregularized problem
% Residuals
r1=L*m1-d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 6: Conjugate Gradient Algorithm
res1 = zeros(M,length(niter)); % residual
res_21 = zeros(length(niter),1); % least-squares residual/objective function
relative_err1 = zeros(length(niter)-1,1);
for ii = 1: length(niter) 
for n = 1:niter(ii)
    gm1 = transpose(L)*r1-lambda*m1;   
    gp1 = transpose(H)*gm1+lambda*p1;
    gm1 = H*gp1;
    gr1 = L*gm1;
    p1 = transpose(gp1)*gp1;
    if n == 1
        sp1 = gp1;
        sm1 = gm1;
        sr1 = gr1;
        beta1 = 0;
    else
        beta1 = p1/p_prep1;
        sp1 = gp1+beta1*sp1;
        sm1 = gm1+beta1*sm1;
        sr1 = gr1+beta1*sr1;
    end
    p_prep1 = p1;
    
    alpha1 = p1/(transpose(sr1)*sr1+lambda*(transpose(sp1)*sp1-transpose(sm1)*sm1));
    p1 = p1-alpha1*sp1;
    m1 = m1-alpha1*sm1;
    r1 = r1-alpha1*sr1;
end
res1(:,ii) = r1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 7: Update the least-squares norm of residuals and relative error
res_21(ii) = transpose(res1(:,ii))*res1(:,ii);
end

% Relative error
for i = 1:length(res_21)-1
    relative_err1(i) = abs(res_21(i+1)-res_21(i))/res_21(i);
end
figure(5)
plot(1:length(m_real), m_real,'r',1:length(m1), m1, '--b')
legend('real model','estimated model')
xlabel('indice of model parameter')
ylabel('model parameter')
prepfig

figure(6)
plot(niter,res_21,'-o')
xlabel('number of iterations')
ylabel('least-squares error')
xlim([min(niter) max(niter)])
prepfig

figure(7)
plot(niter(2:end),relative_err1,'-o')
xlabel('number of iterations')
ylabel('Relative error')
xlim([min(niter) max(niter)])
prepfig


%% Scenario 2 Shaping regularization
p2 = p1;
% Residuals
r2=L*m2-d;

% Step 6: Conjugate Gradient Algorithm
res2 = zeros(M,length(niter)); % residual
res_22 = zeros(length(niter),1); % least-squares residual/objective function
relative_err2 = zeros(length(niter)-1,1);
for ii = 1: length(niter) 
for n = 1:niter(ii)
    gm2 = transpose(L)*r2-lambda*m2; 
    gp2 = transpose(H)*gm2+lambda*p2;
    gm2 = H*gp2;
    gr2 = L*gm2;
    p2 = transpose(gp2)*gp2;
    if n == 1
        sp2 = gp2;
        sm2 = gm2;
        sr2 = gr2;
        beta2 = 0;
    else
        beta2 = p2/p_prep2;
        sp2 = gp2+beta2*sp2;
        sm2 = gm2+beta2*sm2;
        sr2 = gr2+beta2*sr2;
    end
    p_prep2 = p2;
    
    alpha2 = p2/(transpose(sr2)*sr2+lambda*(transpose(sp2)*sp2-transpose(sm2)*sm2));
    p2 = p2-alpha2*sp2;
    m2 = m2-alpha2*sm2;
    r2 = r2-alpha2*sr2;
end
res2(:,ii) = r2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 7: Update the least-squares norm of residuals and relative error
res_22(ii) = transpose(res2(:,ii))*res2(:,ii);
end

% Relative error
for i = 1:length(res_22)-1
    relative_err2(i) = abs(res_22(i+1)-res_22(i))/res_22(i);
end
figure(8)
plot(1:length(m_real), m_real,'r',1:length(m2), m2, '--b')
legend('real model','estimated model')
xlabel('indice of model parameter')
ylabel('model parameter')
prepfig

figure(9)
i = 1:length(niter);
plot(niter,res_22,'-o')
xlabel('number of iterations')
ylabel('least-squares error')
xlim([min(niter) max(niter)])
prepfig

figure(10)
plot(niter(2:end),relative_err2,'-o')
xlabel('number of iterations')
ylabel('Relative error')
xlim([min(niter) max(niter)])
prepfig


%% Scenario 3: Model-reparameterization
d3 = d;
Wd = ones(size(d3));
m3 = zeros(length(m_real),1); 
x3 = m3;
iter = 20;
tol = 1e-10; 
% Wm=(ones(size(m3)));
% x3=myweightcgls(L,d3, Wm,Wd,iter, tol, 1);
for eiter = 1:20
   Wm = abs(x3)+0.1;
   x3 = myweightcgls(L, d3, Wm, Wd, iter, tol, 1);
end

figure(11)
plot(1:length(m_real), m_real,'r',1:length(x3), x3, '--b')
legend('real model','estimated model (model reparameterization)')
xlabel('indice of model parameter')
ylabel('model parameter')
prepfig

%% Scenario 4: Truncated SVD
[U, S, V] = svd(L);
Sp_n = inv(S);
Sp_n(480:end, :) = 0;
%Sp_n(493:end, :) = 0;
Lg_n = V* Sp_n*transpose(U);
m4 = Lg_n * d;

figure(12)
plot(1:length(m_real), m_real,'r',1:length(m4), m4, '--b')
legend('real model','estimated model (truncated SVD)')
xlabel('indice of model parameter')
ylabel('model parameter')
prepfig
end
function [x,rho,eta,nit] = myweightcgls(L,d3, Wm,Wd,k,tol,step)
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
  if (rho(j) < tol) return;end
        
end
return;
end


