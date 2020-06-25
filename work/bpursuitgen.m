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
% Thehe ill conditioned problem to solve is 
% ||d - A x- AH x2||_2^2 + ||x||^1_1 + ||x2||^1_1
%  
% Daniel Trad- UBC - 21-04-2000

clear;
close all;
noisep=5;
freq=40;
dt=0.004;
nt=128;n=nt;
t=0:nt-1;t=t*dt;
niter=4;
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
subplot(311);plot(t,b);title('(a) Signal');ylabel(ylab)
subplot(312);plot(t,e,t,noise);title('(b) Noise');ylabel(ylab)
subplot(313);plot(t,d);title('(d) Data+noise');xlabel(xlab);ylabel(ylab)
end

axis1=[1,length(x)+length(x2),-1.2,1.2];
x_L=zeros(2*nt,1);

tol=1e-10;
eps1=1.5e-5;

AE=[A AH];
de=[d;zeromodel];
xt=[x;x2];
irlstest=0;
karmarkar=0;
nlcg3=1;;

if (irlstest)
  figure(3);
  irls='wpcgnr';
  irls='wtcgls';
  tol=1e-3
  iter_ext=20;
  itercg=nt;
  eps1=1e-1;
  eps2=1e-3;
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'rx'),axis(axis1), title(['IRLS',irls]);hold ...
      on
  tic
  % First pass: predict the most of the noise
  x_irls=IRLS(sparse(AE),d,tol,x0,iter_ext,nt,eps1,eps2,3,irls);
  % Second pass: substarct the noise and compute again
  if (1)
    %tol=1e-5;eps2=1e-1;
    % Apply threshold 
    I=find(abs(x_irls)<0.4);x_irls(I)=0;
    ep=AH*x_irls(length(x)+1:length(xt));dclean=d-ep;
    %keyboard;
    irls='wpcgnr';
    x_irls2=IRLS(sparse(A),dclean,tol,x0(1:nt,1),iter_ext,nt,eps1,eps2,3,irls);
    x_irls=[x_irls2;x_irls(nt+1:2*nt)];
  end
  %%%%%
  
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
  iter_end=25;
  if (noisep==2.5)
    PI_t=1e-5;
    DI_t=1e-5;
    DG_t=1e-5;
    gamma=5e-2;
    delta=1e-1;
    tol=1e-6;
    maxit=70;
  elseif (noisep==5)
    PI_t=1e-5;
    DI_t=1e-5;
    DG_t=1e-5;
    gamma=0.5e-1;
    delta=10e-1;
    tol=1e-5;
    maxit=70;
  end  
    
  figure(4);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), title(['Karmarkar']);
  hold on
  tic
  x_lp=karmarkar_interface(d,AE,iter_end,PI_t,DI_t,DG_t,gamma, ...
			     delta,tol,maxit);
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

if (nlcg3)
  fx='l2l1costparam';
  lambda=1e-0;
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






