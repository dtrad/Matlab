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
nt1=nt;
nt2=nt;
t=0:nt-1;t=t*dt;
niter=4;

% 
% Prior information (Cm^{-1/2}=1/(abs(x)+eps2)+eps1;


w=rickerm(freq,dt);
nw=length(w);

%h=zeros(nw,1);
%h(round(nw/3+1):2*round(nw/3))=1;
h=ones(nw,1);


x1=zeros(nt,1);
x2=zeros(nt1,1);
x3=zeros(nt2,1);

x1(20)=1;
x1(25)=-.9;
x1(30)=-.3;
x1(40)=1;
x1(43)=-0.3;
x1(50)=0.3;
x1(58)=0.7;
x1(63)=-0.1;
x1(80)=0.7;
x1(95)=-0.4;
x1(100)=0.9;

x1=x1(1:nt);

x2(round(nt1/2))=1;
x2(round(nt1/5))=-1;
x2(round(4*nt1/5))=1;

x3(round(nt2/3.2))=-0.5;
x3(round(nt2/2.3))=1;
x3(round(nt2/1.1))=-0.3;


nmodel=length(x1)+length(x2)+length(x3);
zeromodel=zeros(size([x1;x2;x3]));
x0=zeromodel;  

ww=padzeros(w,nt);
hh=padzeros(h,nt1);
% kernel 1
A=circul(ww);
% Kernel 2
AH1=circul(hh);
% Kernel 3
AH2=kernel3(nt2,dt/4);


% data for kernel 1
b=A*x1;
%noise if required.
randn('seed',41997);
%Gaussian noise
e0 =(noise/100)*2*randn(size(b));
%Blocky noise
e1=AH1*x2;
%Harmonic noise
e2=AH2*x3;

d=b+e1+e2;
xlab='time(sec)';
ylab='Amplitude';

if (1)
figure(1)
subplot(321);plot(t(1:nw),w);title('(a) Ricker wavelet');xlabel(xlab);ylabel(ylab)
subplot(322);plot(t,x1);title('(b) Reflectivity');xlabel(xlab);ylabel(ylab)
subplot(323);plot(t(1:nw+40),[zeros(20,1);h;zeros(20,1)]);title('(c) Box car');xlabel(xlab);ylabel(ylab)
subplot(324);plot(t,x2);title('(d) Noise model');xlabel(xlab);ylabel(ylab)
subplot(325);plot(AH2(:,1:50:nt));title('(c) Box car');xlabel(xlab);ylabel(ylab)
subplot(326);plot(t,x3);title('(d) Noise model');xlabel(xlab);ylabel(ylab)


figure(2)
subplot(411);plot(t,b);title('(a) Data');xlabel(xlab);ylabel(ylab)
subplot(412);plot(t,e1);title('(b) Blocky Noise');xlabel(xlab);ylabel(ylab)
subplot(413);plot(t,e2);title('(d) Harmonic noise');xlabel(xlab); ...
    ylabel(ylab)
subplot(414);plot(t,d);title('(d) Data+noise');xlabel(xlab); ...
    ylabel(ylab)

end


% Conjugate gradient solution
% Preconditioned

axis1=[1,length(x1)+length(x2)+length(x3),-1.2,1.2];
x_L=zeros(2*nt,1);

tol=1e-10
iter_ext=1;
iter_ext2=10;
itercg=nt;
eps1=1.5e-5;

%AE=[[A AH];eps1*eye(2*nt)];
AE=[A AH1 AH2];
de=[d;zeromodel];

%
irlstest=1;
karmarkar=0;
nlcg1=0;
nlcg2=0;

if (irlstest)
  figure(3);
  irls='wpcgnr';
  eps1=1e-1;
  eps2=1e-3;
  subplot(1,1,1)
  plot(1:nmodel,[x1;x2;x3],'rx'),axis(axis1), title(['IRLS',irls]);hold ...
      on
  x_irls=IRLS(sparse(AE),d,tol,x0,10,nt,eps1,eps2,3,irls);
  plot(1:length(x_irls),x_irls,'r'),axis(axis1), title(['IRLS',irls]);
  hold off
  figure(10)
  subplot(311)
  plot(1:nmodel,x_irls,1:nmodel,[x1;x2;x3],'x'),     
  axis(axis1),xlabel('L = I'), title(['IRLS',irls])
end

if (karmarkar)
  % Karmarkar interior point method
  iter_end=100;
  PI_t=1e-3;
  DI_t=1e-3;
  DG_t=1e-3;
  gamma=1;
  delta=1e-3;
  figure(4);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), title(['Karmarkar']);
  hold on
  x_cgnl=karmarkar_interface(d,AE,iter_end,PI_t,DI_t,DG_t,gamma,delta);
  hold off;
  figure(10)
  subplot(312)
  plot(1:length(x_cgnl),x_cgnl,1:length([x;x2]),[x;x2],'x'),axis(axis1), title(['karmarkar']);

  hold off
end

if (nlcg1)
  fx='l2l1cost';
  lambda=1e-1;
  itercg=3*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  [x_cgnl1,delnew,niter]=nlcgproject4(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,0);
  figure(10);
  subplot(3,1,3)
  plot(1:nmodel,x_cgnl1,1:nmodel,[x;x2],'x'),axis(axis1), title(['nlcg']);
  
end

if (nlcg2)
  fx='l2l1cost'
  lambda=1e-2;
  itercg=3*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  [x_cgnl5,fy,normdeltay,normgrad]=nlcgproject5(fx,x0,d,sparse(AE),lambda,itercg,tol); % NL CG
  %[x_cgnl,delnew,niter]=nlcgproject4(d1fx,d2fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,0);
  figure(10);
  subplot(3,1,3)
  plot(1:nmodel,x_cgnl5,1:nmodel,[x;x2],'x'),axis(axis1), title(['nlcg']);
  
end

return;






