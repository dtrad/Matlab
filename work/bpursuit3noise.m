% program to test algorithms for sparse inversion 
% The signal is  a sparse reflectivity series convolved with a Ricker
% wavelet (deconvolution with known wavelet). The data contains
% noise in the form of coherent and uncoherent (random).
% The inversion is performed with the purpose of modeling signal and noise
% at the same time. The hybrid model can then used to predict the
% signal or noise. 
% This wavelet is included in the kernel A;
% The coherent noise is included in the kernel AH;
% Both are  a circulant matrices, 
% The ill conditioned problem to solve is 
% ||d - A x- AH x2||_2^2 + ||x||^1_1 + ||x2||^1_1
%  
% Daniel Trad- UBC - 10-03-2001

clear;
close all;
noisep=5;
freq=40;
dt=0.004;
nt=128;n=nt;
t=0:nt-1;t=t*dt;
niter=4;

% Prior information (Cm^{-1/2}=1/(abs(x)+eps2)+eps1;
w=rickerm(freq,dt);
nw=length(w);
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
if (0) AH=kernel3(nt,dt/4);end
% data
b=A*x;
%noise
randn('seed',41997);
%Gaussian noise
noise =(noisep/100)*2*randn(size(b));
%Blocky noise
e=AH*x2;

d=b+e+noise;
xlab='time(sec)';
ylab='Amplitude';

if (1)
figure(1)
subplot(221);plot(t(1:nw),w);title('(a) Ricker wavelet');ylabel(ylab)
subplot(222);plot(t,x);title('(b) Reflectivity');ylabel(ylab)
subplot(223);plot(t(1:nw+40),[zeros(20,1);h;zeros(20,1)]);title('(c) Box car');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,x2);title('(d) Noise model');xlabel(xlab);ylabel(ylab)
figure(2)
axis1=[0 2*nt 1.5 1.5];
axis2=[0 nt min(d) max(d)];
axis3=[0 t(nt) min(d) max(d)];
subplot(311);plot(t,b);title('(a) Signal');ylabel(ylab);axis(axis3);
subplot(312);plot(t,e,t,noise);title('(b) Noise');ylabel(ylab);axis(axis3);
subplot(313);plot(t,d);title('(d) Data+noise');xlabel(xlab); 
ylabel(ylab);axis(axis3);
end

% Conjugate gradient solution
% Preconditioned
axis1=[1,length(x)+length(x2),-1.2,1.2];
tol=1e-10;
eps1=1.5e-5;

AE=[A AH];
de=[d;zeromodel];
xt=[x;x2];

irlstest=0;  
karmarkar=0;
nlcg1=0;
nlcg2=0;
nlcg3=1;



if (irlstest)
  figure(3);

  irls='wpcgnr';
  %irls='wtcgls';
  method=['Log Barrier',irls]
  tol=1e-14
  iter_ext=5;
  itercg=2*nt;
  if (noisep==0) eps1=5e-2; eps2=1e-3;
  elseif (noisep==2.5) eps1=30e-2; eps2=1e-2;
  elseif (noisep==5) eps1=30e-2; eps2=2e-2;
  end 
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'rx'),axis(axis1), title(['IRLS', irls]);
  %hold on
  tic
  % First pass: predict the most of the noise
  %x_sol=IRLS(sparse(AE),d,tol,x0,iter_ext,itercg,eps1,eps2,3,irls);
  x_sol=IRLSext(sparse(AE),d,tol,x0,iter_ext,itercg,eps1,eps2,3,irls);
  % Second pass: substarct the noise and compute again
  if (0)
    % Apply threshold 
    % I=find(abs(x_irls)<0.4);x_irls(I)=0;    
    ep=AH*x_sol(length(x)+1:length(xt));dclean=d-ep;
    x_sol2=IRLSext(sparse(A),dclean,tol,x0(1:nt,1),iter_ext,itercg/2,eps1,eps2,3,irls);
    %x_sol2=IRLS(sparse(A),dclean,tol,x0(1:nt,1),iter_ext,nt,eps1,eps2,3,irls);
    x_sol=[x_sol2;x_sol(nt+1:2*nt)];
  end
  %%%%%
end

if (karmarkar)
  % Karmarkar interior point method
  method=['Log Barrier']
  iter_end=25;
  if (noisep==0)
    iter_end=15;
    PI_t=1e-5;
    DI_t=1e-5;
    DG_t=1e-5;
    gamma=1e-7;
    delta=1e-7;
    tol=1e-10;
    maxit=200;
  elseif (noisep==2.5)
    PI_t=1e-5;
    DI_t=1e-5;
    DG_t=1e-5;
    gamma=5e-2;
    delta=1e-1;
    tol=1e-6;
    maxit=70;
  elseif (noisep==5)
    iter_end=15;
    PI_t=1e-6;
    DI_t=1e-6;
    DG_t=1e-6;
    gamma=0.001e-1;
    delta=6e-1;
    tol=1e-6;
    maxit=150;
  end  
    
  figure(4);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1),title(['Log Barrier ']);
  hold on
  tic
  x_sol=karmarkar_interface(d,AE,iter_end,PI_t,DI_t,DG_t,gamma, ...
			     delta,tol,maxit);
end

if (nlcg1)
  method=['Nonlinear CG']
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
  [x_sol,delnew,niter]=nlcgproject4(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,0,xt);
  
end

if (nlcg2)
  method=['Nonlinear CG']
  fx='l2l1cost'
  lambda=1e-2;
  itercg=3*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  [x_sol,fy,normdeltay,normgrad]=nlcgproject5(fx,x0,d,sparse(AE),lambda,itercg,tol); % NL CG
  
end
if (nlcg3)
  method=['Nonlinear CG']
  fx='l2l1costparam';
  if (noisep==0) lambda=1e-0;tol=1e-5;itercg=100
  elseif (noisep==2.5) lambda=5e-2;tol=1e-2;itercg=200
  elseif (noisep==5) lambda=50e-2;tol=4e-2;itercg=200
  end
  %lambda=1e-0;
  
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  tic
  [x_sol,delnew,niter]=nlcgproject6(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,1,xt);
    
end

  
text1=sprintf('Noise=%4.1f %% , required time %5.2f s\n',noisep, toc);
dp=AE*x_sol;
bp=A*x_sol(1:length(x));
ep=AH*x_sol(length(x)+1:length(xt));

hold off  
plot(1:length(x_sol),x_sol,'r'),axis(axis1), 

figure(10)
subplot(311)
plot(1:nmodel,x_sol,1:nmodel,[x;x2],'.'),     
axis(axis1),title([method,' : ',text1])
subplot(312)
plot(1:length(b),b,'.',1:length(b),bp),axis(axis2);  
title('signal and predicted signal')    
subplot(313)
plot(1:length(e),e,'.',1:length(e),ep),axis(axis2);  
title('noise and predicted noise')    


return;






