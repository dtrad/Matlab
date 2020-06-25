% program to test hyperparameter for sparse deconvolution 
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
% Daniel Trad- UBC - 21-04-2000

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


figure(1)
subplot(221);plot(t(1:nw),w);title('(a) Ricker wavelet');xlabel(xlab);ylabel(ylab)
subplot(222);plot(t,x);title('(b) Reflectivity');xlabel(xlab);ylabel(ylab)
subplot(223);plot(t,b);title('(c) Data');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,d);title('(d) Data+noise');xlabel(xlab);ylabel(ylab)


% Conjugate gradient solution
% Preconditioned

axis1=[1,n,-1.2,1.2];
x_L=zeros(nt,1);
tol=1e-2
iter_ext=10;
itercg=13;

[x_I,rho,niter]=wpcgnr2(A,d,eye(nt),eye(nt),tol,zeros(nt,1), ...
			itercg);

figure(2)
subplot(3,1,1);
plot(1:n,x_I,1:n,x,'x'),axis(axis1),xlabel('L = I'), title(['zero' ...
		    ' order regularization'])

x_L=x_I;
eps1=1e-3;
ifig=1;

for iter=1:iter_ext
  %L=diag(1./(max(abs(x_L),eps1)));
  L=100*diag((abs(x_L))+0.01);
  [x_L,rho,niter,J(iter)]=wpcgnr2(A,d,L,eye(nt),tol,zeros(nt,1),itercg);
  display('iter CG    iter #   J');
  [niter iter J(iter)]

  if (mod(iter,5)==0)&(iter>1)
    ifig=ifig+1;
    figure(2)
    mytext=sprintf('iter=%d\n',iter);
    subplot(3,1,ifig);
    plot(1:n,x_L,1:n,x,'x'),axis(axis1),xlabel('L \neq I'),ylabel(mytext)
    axis(axis1);
  end
  %if (mod(iter,2) ==0 ) itercg=itercg-1; end
  %if (iter==iter_ext-1) itercg=7;end
end;
xsol=x_L;

figure(3),
subplot(211);plot(x);title('(a) Initial model');
subplot(212);plot(xsol);
mytext=sprintf('(b) Weighted solution %d iterations\n',iter);
title(mytext);

figure(4)
plot(J)








