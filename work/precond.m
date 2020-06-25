function [x_g,x_g2]=precond(N,beta)
% Regulariztion and standard form
% Define A and Wm and look at the svd for different Hessian forms
if (nargin < 1) N=10;M=N;end
if (nargin < 2) beta=1;end

A=randn(N,N)+eye(N);
x=randn(N,1);
b=A*x;
e=ones(N,1);
Wm=diag(ones(N-1,1),1) + diag(ones(N-1,1),-1) - 2*diag(ones(N,1),0);
%W=diag([e -2*e e],0:1);
AA=A'*A;
H1=AA+beta*Wm'*Wm;
Wmi=inv(Wm);
H2=Wmi'*AA*Wmi+beta*eye(N);
[A_s,b_s,Wm_p,K,M] = std_form(A,Wm,b);
H3=A_s'*A_s + beta * eye(N);  
figure,
subplot(221);semilogy(sqrt(svd(AA))); V=axis;
subplot(222);semilogy(sqrt(svd(H1))); axis(V);
subplot(223);semilogy(sqrt(svd(H2))); axis(V);
subplot(224);semilogy(sqrt(svd(H3))); axis(V);

% Stage 3
% Compare iterations
display('General form iter=')
[x_g,rho,eta] = wtcgls(A,Wm,b,N,0);
display('Standard Form')
[x_s,rho,eta] = wtcgls(A_s,eye(N),b_s,N,0);
x_g2 = gen_form(Wm_p,x_s,zeros(size(x_s)));
figure,
i=1:100;
subplot(211),plot(i,x,i,x_g);
subplot(212),plot(i,x,i,x_g2); 











