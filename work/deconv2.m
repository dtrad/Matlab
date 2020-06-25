% 
% 'program to test sparse deconvolution 
% The wavelet is deconvolved with CG
% A given Ricker wavelet is used for the example.
% This wavelet is included in the kernel A;
% which is a circulant matrix, to perform convolution on the
% unknown model (here is given) to produce the data.
% This means that the ill conditioned problem to solve is 
% ||d-Ax||^2 + ||Lx||^2
% where L is our prior
% For this example I used the L1 prior
% L=1/abs(x);
% to get sparseness.
%  
% Based in Hanke s regularization toolbox
% Daniel Trad- UBC - 13-07-99'
clear;
close all;
noise=5;
freq=40;
dt=0.004;
nt=128;n=nt;
t=0:nt-1;t=t*dt;
niter=4;

% 
% Prior information (Cm^{-1/2}=1/(abs(x)+eps2)+eps1;

eps1=1e-5;  
eps2=1e-5;

w=rickerm(freq,dt);
nw=length(w);
ww=[w((nw+1)/2:nw)];
x=zeros(nt,1);
x(20)=1;
x(25)=-.9;
x(30)=-.3;
x(40)=1;
x(43)=-0.3;
x(50)=0.3;
x(58)=0.7;
x(63)=-0.1;
x(80)=0.7;
x(95)=-0.4;
x(100)=0.9;

ww=padzeros(ww,nt);

% kernel
A=toeplitz(ww);
% data
b=A*x;
%noise
randn('seed',41997);
e =(noise/100)*2*randn(size(b));
d=b+e;
xlab='time(sec)';
ylab='Amplitude';


figure(1)
subplot(221);plot(t(1:nw),w);title('Ricker wavelet');xlabel(xlab);ylabel(ylab)
subplot(222);plot(t,x);title('Reflectivity');xlabel(xlab);ylabel(ylab)
subplot(223);plot(t,b);title('Reflectivity+noise');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,d);title('solution');xlabel(xlab);ylabel(ylab)

% Conjugate gradient solution

b=d;  
x_0=zeros(nt,1);
tol=1e-2
itercg=20;
beta=1e-3;
eps1=1e-3;
L=diag(1./(max(abs(x),eps1)));
N=n;
display('=====> General form iter')
% True regularization
[x_g1,rho_g1,eta_g1,iter_g1] = wtcgls([A;sqrt(beta)*L],eye(n),[b; ...
		    zeros(size(b))],n,tol);

% Solution in the subspace
[x_g2,rho_g2,eta_g2,iter_g2] = wtcgls(A,sqrt(beta)*L,b,n,tol);

display('=====>Standard Form')

% transformation to standard form is easy here because Wm is square
% and full range.
L_p = inv(L); A_s = A/L; b_s = b;
% For more complicated cases use std_form as follows
% [A_s,b_s,Wm_p,K,M] = std_form(A,Wm,b);

[x_s,rho_s1,eta_s1,iter_s1] = wtcgls([A_s;sqrt(beta)*eye(n)], ...
				     eye(N),[b_s;zeros(size(b_s))],n,tol);
x_s1 = L_p * x_s;
%x_s1 = gen_form(Wm_p,x_s,zeros(size(x_s)));
% Solution in the subspace
[x_s,rho_s2,eta_s2,iter_s2] = wtcgls(A_s,eye(n),b_s,n,tol);
% For this simple case use
x_s2 = L_p * x_s;
% For more complicated cases use
% x_s2 = gen_form(Wm_p,x_s,zeros(size(x_s)));
i=1:length(x);
figure,
subplot(221),plot(i,x,'+',i,x_g1);title('General Form, regularization');
subplot(222),plot(i,x,'+',i,x_g2);title('General Form, subspace'); 
subplot(223),plot(i,x,'+',i,x_s1);title('Standard Form, regularization');
subplot(224),plot(i,x,'+',i,x_s2);title('Standard Form, subspace');
figure,
subplot(221),semilogy(rho_g1);title('Residuals, General Form, regularization');
subplot(222),semilogy(rho_g2);title('Residuals Form, subspace');
subplot(223),semilogy(rho_s1);title('Residuals Standard Form, regularization');
subplot(224),semilogy(rho_s2);title('Residuals Standard Form, subspace'); 

























