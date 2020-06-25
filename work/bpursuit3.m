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


w=rickerm(freq,dt);
nw=length(w);

%h=zeros(nw,1);
%h(round(nw/3+1):2*round(nw/3))=1;
h=ones(nw,1);


x=zeros(nt,1);
x2=zeros(nt,1);

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

x2(round(nt/2))=1;
x2(round(nt/5))=-1;
x2(round(4*nt/5))=1;

nmodel=length(x)+length(x2);
zeromodel=zeros(size([x;x2]));
x0=zeromodel;  

ww=padzeros(w,nt);
hh=padzeros(h,nt);
% kernel
A=circul(ww);
AH=circul(hh);
% data
b=A*x;
%noise
randn('seed',41997);
%Gaussian noise
e =(noise/100)*2*randn(size(b));
%Blocky noise
e=AH*x2;

d=b+e;
xlab='time(sec)';
ylab='Amplitude';

if (1)
figure(1)
subplot(221);plot(t(1:nw),w);title('(a) Ricker wavelet');xlabel(xlab);ylabel(ylab)
subplot(222);plot(t,x);title('(b) Reflectivity');xlabel(xlab);ylabel(ylab)
subplot(223);plot(t(1:nw+40),[zeros(20,1);h;zeros(20,1)]);title('(c) Box car');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,x2);title('(d) Noise model');xlabel(xlab);ylabel(ylab)
figure(2)
subplot(311);plot(t,b);title('(a) Data');xlabel(xlab);ylabel(ylab)
subplot(312);plot(t,e);title('(b) Noise');xlabel(xlab);ylabel(ylab)
subplot(313);plot(t,d);title('(d) Data+noise');xlabel(xlab);ylabel(ylab)
end

% Conjugate gradient solution
% Preconditioned

axis1=[1,length(x)+length(x2),-1.2,1.2];
x_L=zeros(2*nt,1);

tol=1e-10;
eps1=1.5e-5;

%AE=[[A AH];eps1*eye(2*nt)];
AE=[A AH];
de=[d;zeromodel];
xt=[x;x2];
%
irlstest=0;
karmarkar=0;
nlcg1=0;
nlcg2=0;
nlcg3=1;;

if (irlstest)
  figure(3);
  irls='wpcgnr';
  %irls='wtcgls';
  tol=1e-10
  iter_ext=10;
  itercg=nt;
  eps1=1e-1;
  eps2=1e-7;
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'rx'),axis(axis1), title(['IRLS',irls]);hold ...
      on
  tic
  x_irls=IRLS(sparse(AE),d,tol,x0,iter_ext,nt,eps1,eps2,3,irls);
  text1=sprintf('required time %5.2f s\n',toc);
  dp=AE*x_irls;
  bp=A*x_irls(1:length(x));
  ep=AH*x_irls(length(x)+1:length(xt));
  
  plot(1:length(x_irls),x_irls,'r'),axis(axis1), title(['IRLS',irls]);
  hold off
  figure(10)
  subplot(311)
  plot(1:nmodel,x_irls,1:nmodel,[x;x2],'x'),     
  axis(axis1),title(['IRLS ',irls,' : ',text1])
  subplot(312)
  plot(1:length(b),b,1:length(b),bp,'.'),  
  title('signal and predicted signal')    
  subplot(313)
  plot(1:length(e),e,1:length(e),ep,'.'), 
  title('noise and predicted noise')    
end

if (karmarkar)
  % Karmarkar interior point method
  iter_end=15;

  PI_t=1e-5;
  DI_t=1e-5;
  DG_t=1e-5;
  gamma=1e-7;
  delta=1e-7;
  
  figure(4);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), title(['Karmarkar']);
  hold on
  tic
  x_lp=karmarkar_interface(d,AE,iter_end,PI_t,DI_t,DG_t,gamma, ...
			     delta);
  text1=sprintf('required time %5.2f s\n',toc);
  dp=AE*x_lp;
  bp=A*x_lp(1:length(x));
  ep=AH*x_lp(length(x)+1:length(xt));

  hold off;
  figure(10)
  subplot(311)
  plot(1:length(x_lp),x_lp,1:length([x;x2]),[x;x2],'x'),axis(axis1), title(['karmarkar']);
  axis(axis1),title(['Karmarkar',text1])
  subplot(312)
  plot(1:length(b),b,1:length(b),bp,'.'),  
  title('signal and predicted signal')    
  subplot(313)
  plot(1:length(e),e,1:length(e),ep,'.'), 
  title('noise and predicted noise')    
  hold off
end

if (nlcg1)
  fx='l2l1cost';
  lambda=1e-1;
  tol=1e-5;
  itercg=1000;%*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  tic
  [x_cgnl1,delnew,niter]=nlcgproject4(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,0,xt);
  text1=sprintf('required time %5.2f s\n',toc);
  dp=AE*x_cgnl1;
  bp=A*x_cgnl1(1:length(x));
  ep=AH*x_cgnl1(length(x)+1:length(xt));

  figure(10);
  subplot(311)
  plot(1:nmodel,x_cgnl1,1:nmodel,[x;x2],'x'),axis(axis1), title(['nlcg',text1]);
  subplot(312)
  plot(1:length(b),b,1:length(b),bp,'.'),  
  title('signal and predicted signal')    
  subplot(313)
  plot(1:length(e),e,1:length(e),ep,'.'), 
  title('noise and predicted noise')    
  hold off
  
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
if (nlcg3)
  fx='l2l1costparam';
  lambda=1e-1;
  tol=1e-5;
  itercg=500;%*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  tic
  [x_cgnl1,delnew,niter]=nlcgproject6(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,0,xt);
  text1=sprintf('required time %5.2f s\n',toc);
  dp=AE*x_cgnl1;
  bp=A*x_cgnl1(1:length(x));
  ep=AH*x_cgnl1(length(x)+1:length(xt));

  figure(10);
  subplot(311)
  plot(1:nmodel,x_cgnl1,1:nmodel,[x;x2],'x'),axis(axis1), title(['nlcg',text1]);
  subplot(312)
  plot(1:length(b),b,1:length(b),bp,'.'),  
  title('signal and predicted signal')    
  subplot(313)
  plot(1:length(e),e,1:length(e),ep,'.'), 
  title('noise and predicted noise')    
  hold off
  
end
return;






