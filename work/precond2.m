function [A,b,x]=precond2(N,beta,tol,stage)
% Input 
%       N size
%       beta hyperparameter
%       tol  tolerance for (res'*res)       
%       stage 3,4  Algorithm to test
% Regularization and standard form
% Define A and Wm and look at the svd for different Hessian forms
% Daniel Trad ----- UBC 

if (nargin < 1) N=10;M=N;end
if (nargin < 2) beta=1;end
if (nargin < 3) tol=1e-7;end
if (nargin < 4) stage=3;end
 
i=1:N;
x=linspace(0,1,N);
for ii=1:N 
A(ii,:)=exp(-ii*x/2).*cos(ii*x/2);
end
figure,
subplot(211),plot(A');
%A=randn(N,N)+eye(N);
dx=round(N/4);
x=[zeros(dx,1);ones(dx,1);0.5*ones(dx,1);zeros(dx,1)];
x=x(1:N);
subplot(212),plot(x);
 
b=A*x;
e=ones(N,1);
Wm=diag(ones(N-1,1),1) + diag(ones(N-1,1),-1) - 2*diag(ones(N,1),0);
%W=diag([e -2*e e],0:1);
AA=A'*A;
H1=AA+beta*Wm'*Wm;
Wmi=inv(Wm);
H2=Wmi'*AA*Wmi+beta*eye(N);
% transformation to standard form is easy here because Wm is square
% and full range.
Wm_p = inv(Wm); A_s = A/Wm; b_s = b;
% For more complicated cases use std_form as follows
% [A_s,b_s,Wm_p,K,M] = std_form(A,Wm,b);

H3=A_s'*A_s + beta * eye(N);  
figure,
subplot(221);semilogy(sqrt(svd(AA))); V=axis;
subplot(222);semilogy(sqrt(svd(H1))); axis(V);
subplot(223);semilogy(sqrt(svd(H2))); axis(V);
subplot(224);semilogy(sqrt(svd(H3))); axis(V);

if (stage==3)
% Stage 3
% Compare iterations
display('General form iter=')
% True regularization
[x_g1,rho_g1,eta_g1,iter_g1] = wtcgls([A;sqrt(beta)*Wm],eye(N),[b; ...
		    zeros(size(b))],N,tol);

% Solution in the subspace
[x_g2,rho_g2,eta_g2,iter_g2] = wtcgls(A,sqrt(beta)*Wm,b,N,tol);

display('Standard Form')
[x_s,rho_s1,eta_s1,iter_s1] = wtcgls([A_s;sqrt(beta)*eye(N)], ...
				     eye(N),[b_s;zeros(size(b_s))],N,tol);
x_s1 = Wm_p * x_s;
%x_s1 = gen_form(Wm_p,x_s,zeros(size(x_s)));
% Solution in the subspace
[x_s,rho_s2,eta_s2,iter_s2] = wtcgls(A_s,eye(N),b_s,N,tol);
% For this simple case use
x_s2 = Wm_p * x_s;
% For more complicated cases use
% x_s2 = gen_form(Wm_p,x_s,zeros(size(x_s)));

figure,

subplot(221),plot(i,x,i,x_g1);
subplot(222),plot(i,x,i,x_g2); 
subplot(223),plot(i,x,i,x_s1);
subplot(224),plot(i,x,i,x_s2);
figure,
subplot(221),semilogy(rho_g1);
subplot(222),semilogy(rho_g2);
subplot(223),semilogy(rho_s1);
subplot(224),semilogy(rho_s2); 

%[iter_g1 iter_g2 iter_s1 iter_s2];

elseif (stage==4)

% Stage 4 
% Comparison of Methods with Precondtioning
% Here we try some preconditioners

%M=diag(sqrt(diag(Wmi'*AA*Wmi)));
M=eye(N);
%M=diag(sqrt(diag(AA)));
%M=diag(sqrt(diag(Wmi'*AA*Wmi+beta*eye(N))));
%M=diag(diag(H1));
MS=diag(diag(A'*A));

% True regularization

% WTCGLS
[x_g3,rho_g3,eta_g3,iter_g3] = wtcgls([A;sqrt(beta)*Wm],M,[b;zeros(size(b))],N,tol);

% LPCGNR
[x_g4,rho_g4,iter_g4]=lpcgnr([A;sqrt(beta)*Wm],[b;zeros(size(b))],M,tol);

% PCGLS Hanke Toolbox
%[x_s,rho,eta,F] = pcgls(A_s,eye(N),[1;zeros(N-1,1)],b_s,N,0)
%x_s3 = gen_form(Wm_p,x_s,zeros(size(x_s)));

% LPCGNR: Standard form 
[x_s,rho_s3,iter_s3]=lpcgnr([A_s;sqrt(beta)*eye(N)],[b_s;zeros(size(b))],M,tol);
x_s3 = gen_form(Wm_p,x_s,zeros(size(x_s)));

% Applying CG to A'A x =A'b with H=A'A and b=A'b
AE=[A;sqrt(beta)*Wm];
H1=AE'*AE;
[x_g5,rho_g5,iter_g5]=pcgsaad(A'*A,A'*b,MS,tol);

%[x_g5,rho_g5,iter_g5]=pcgsaad(H1,AE'*[b;zeros(size(b))],eye(N),tol);

%AE=[A*Wmi;sqrt(beta)*eye(N)];
%H2=AE'*AE;
%H2=Wmi'*AA*Wmi+beta*eye(N);
% PCG in standard form
%[x_s,rho_s4,iter_s4]=pcgsaad(H2,AE'*[b;zeros(size(b))],eye(N),tol);
%x_s4 = gen_form(Wm_p,x_s,zeros(size(x_s)));

figure,
subplot(221),plot(i,x,i,x_g3);V=axis;
subplot(222),plot(i,x,i,x_g4);axis(V);
subplot(223),plot(i,x,i,x_s3);axis(V);
subplot(224),plot(i,x,i,x_g5);axis(V);

figure
subplot(221),semilogy(rho_g3);V=axis;
subplot(222),semilogy(rho_g4);axis(V);
subplot(223),semilogy(rho_s3);axis(V);
subplot(224),semilogy(rho_g5);axis(V);


[iter_g3 iter_g4 iter_s3 iter_g5]
end
display('Here---------------------------\n');





