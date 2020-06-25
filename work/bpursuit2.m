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

nm=length(x)+length(x2);


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

tol=1e-10
iter_ext=1;
iter_ext2=10;
itercg=nt;
eps1=1.5e-5;

%AE=[[A AH];eps1*eye(2*nt)];
AE=[A AH];
de=[d;zeros(size(x));zeros(size(x2))];

[x_cgnl,delnew,niter]=nlcgproject(AE,d,tol,zeros(2*nt,1),eye(2*nt),1e-3,3*length(x_L));
figure(3);
subplot(4,1,2)
plot(1:length(x_cgnl),x_cgnl,1:length([x;x2]),[x;x2],'x'),axis(axis1), title(['nlcg']);

return;


[x_I,rho,niter]=wpcgnr2(AE,de,eye(2*nt),tol,zeros(2*nt,1), ...
			itercg);

figure(3)
subplot(4,1,1);
plot(1:length(x_I),x_I,1:length([x;x2]),[x;x2],'x'),axis(axis1),xlabel('L = I'), title(['zero' ...
		    ' order regularization'])

x_L=x_I;
eps2=1e-2;
ifig=1;

for iter=1:iter_ext
  %L=diag(1./(max(abs(x_L),eps2)));
  %L=100*diag((abs(x_L))+0.01);
  %L=diag(1./sqrt(x_L.^2+eps2.^2)); 
  L=diag(1./sqrt([x;x2].^2+eps2.^2)); 
  %L=diag(1./sqrt([x;zeros(size(x2))].^2+eps2.^2)); 
  AE=[[A AH];eps1*L];
  [x_L,rho,niter,Jd(iter)]=wpcgnr2(AE,de,eye(2*nt),tol,zeros(2*nt,1),itercg);
  Jm(iter)=(x_L'*((eps1*L)*x_L));
  J(iter)=Jd(iter)+Jm(iter);
  display('iter CG    iter #   J   Jd    Jm'  );
  [niter iter J(iter) Jd(iter) Jm(iter)]

  if ((mod(iter,5)==0)&(iter>1))|(iter==1)

    ifig=ifig+1;
    figure(3)
    mytext=sprintf('iter=%d\n',iter);
    subplot(4,1,ifig);
    plot(1:length(x_L),x_L,1:length([x;x2]),[x;x2],'x'),axis(axis1),xlabel('L \neq I'),ylabel(mytext)
    axis(axis1);
  end
  %if (mod(iter,2) ==0 ) itercg=itercg-1; end
  %if (iter==iter_ext-1) itercg=7;end
end;
x_L=x_I;
for iter=1:iter_ext2
  %L=diag(1./(max(abs(x_L),eps2)));
  %L=100*diag((abs(x_L))+0.01);
  %L=diag(1./sqrt(x_L.^2+eps2.^2)); 
  %L=diag(1./sqrt([x;x2].^2+eps2.^2)); 
  L=diag(1./sqrt([x_L(1:length(x));x2].^2+eps2.^2)); 
  
  AE=[[A AH];eps1*L];
  [x_L,rho,niter,Jd(iter)]=wpcgnr2(AE,de,eye(2*nt),tol,zeros(2*nt,1),itercg);
  Jm(iter)=(x_L'*((eps1*L)*x_L));
  J(iter)=Jd(iter)+Jm(iter);
  display('iter CG    iter #   J   Jd    Jm'  );
  [niter iter J(iter) Jd(iter) Jm(iter)]

  if ((mod(iter,iter_ext2)==0)&(iter>1))|(iter==iter_ext2)

    ifig=ifig+1;
    figure(3)
    mytext=sprintf('iter=%d\n',iter);
    subplot(4,1,ifig);
    plot(1:length(x_L),x_L,1:length([x;x2]),[x;x2],'x'),axis(axis1),xlabel('L \neq I'),ylabel(mytext)
    axis(axis1);
  end
  %if (mod(iter,2) ==0 ) itercg=itercg-1; end
  %if (iter==iter_ext-1) itercg=7;end
end;
x_L=x_I;
for iter=1:iter_ext2
  %L=diag(1./(max(abs(x_L),eps2)));
  %L=100*diag((abs(x_L))+0.01);
  L=diag(1./sqrt(x_L.^2+eps2.^2)); 
  %L=diag(1./sqrt([x;x2].^2+eps2.^2)); 
  %L=diag(1./sqrt([x_L(1:length(x));x2].^2+eps2.^2)); 
  
  AE=[[A AH];eps1*L];
  [x_L,rho,niter,Jd(iter)]=wpcgnr2(AE,de,eye(2*nt),tol,zeros(2*nt,1),itercg);
  Jm(iter)=(x_L'*((eps1*L)*x_L));
  J(iter)=Jd(iter)+Jm(iter);
  display('iter CG    iter #   J   Jd    Jm'  );
  [niter iter J(iter) Jd(iter) Jm(iter)]

  if ((mod(iter,iter_ext2)==0)&(iter>1))|(iter==iter_ext2)

    ifig=ifig+1;
    figure(3)
    mytext=sprintf('iter=%d\n',iter);
    subplot(4,1,ifig);
    plot(1:length(x_L),x_L,1:length([x;x2]),[x;x2],'x'),axis(axis1),xlabel('L \neq I'),ylabel(mytext)
    axis(axis1);
  end
  %if (mod(iter,2) ==0 ) itercg=itercg-1; end
  %if (iter==iter_ext-1) itercg=7;end
end;


xsol=x_L;

figure(4),
subplot(211);plot([x;x2]);title('(a) Initial model');
axis([1 length(xsol) min(xsol) max(xsol)]);
subplot(212);plot(xsol);
mytext=sprintf('(b) Weighted solution %d iterations\n',iter);
title(mytext);
axis([1 length(xsol) min(xsol) max(xsol)]);

figure(5)
plot(J);title('Cost function');xlabel('# iterations');ylabel('J')     







