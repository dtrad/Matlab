% program to test sparse deconvolution 
% The wavelet is deconvolved with nonlinear CG
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
% Daniel Trad- UBC - 21-04-2000
clear;
close all;
noise=5;
freq=40;
dt=0.004;
nt=128;n=nt;
t=0:nt-1;t=t*dt;
niter=4;

% Prior information (Cm^{-1/2}=1/(abs(x)+eps2)+eps1;

eps1=1e-5;  
eps2=1e-5;

w=rickerm(freq,dt);
nw=length(w);

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

ww=padzeros(w,nt);

% kernel
A=circul(ww);
% data
b=A*x;
%noise
randn('seed',41997);
e =(noise/100)*2*randn(size(b));
d=b+e;
xlab='time(sec)';
ylab='Amplitude';

if (0)
figure(1)
subplot(221);plot(t(1:nw),w);title('Ricker wavelet');xlabel(xlab);ylabel(ylab)
subplot(222);plot(t,x);title('Reflectivity');xlabel(xlab);ylabel(ylab)
subplot(223);plot(t,b);title('Reflectivity+noise');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,d);title('solution');xlabel(xlab);ylabel(ylab)
end
% Conjugate gradient solution

%b=d;  
x_0=zeros(nt,1);
tol=1e-4;
itercg=10*nt;
beta=1e-3;
eps1=1e-3;
%L=diag(1./(max(abs(x),1e-3)));
L=eye(length(x));
N=n;
display('=====> Nonlinear CG')
lambda=0.1;
[x_NLCG,delnew,niter]=nlcgNRFR(A,d,tol,x_0,L,lambda,itercg);
delnew
niter

% True regularization
i=1:length(x);

figure,
plot(i,x,'+',i,x_NLCG);title('Non Linear CG');


























