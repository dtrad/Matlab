% program to test sparse deconvolution 
% The wavelet is deconvolved with Tikhonov s regularization,
% and truncated SVD. The best of these is used as prior information
% for TGSVD. 
% Two stopping criterium are used: L-curve and GCV
% A given Ricker wavelet is used for the example.
% This wavelet is included in the kernel A;
% which is a circulant matrix, to perform convolution on the
% unknown model (here is given) to produce the data.
% This means that the ill conditioned problem to solve is 
% Ax=d
% TSVD solves the problem ||d-Ax||^2 + ||x||^2
% GTSD solves the same problem with prior information
% ||d-Ax||^2 + ||Lx||^2
% where L is our prior
% For this example I used the L1 prior
% L=1/abs(x+eps)+eps;
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
subplot(223);plot(t,b);title('(c) data');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,d);title('(d) data+noise');xlabel(xlab);ylabel(ylab)


% Let us studt first the kernel A and the data

[U,s,V] = csvd(A);
display(sprintf('Cond number for A=%f',cond(A)));
W=fft(padzeros(w,nt));
WW=W.*conj(W);
f=(-nt/2):(nt/2-1);f=f/(nt*dt);
figure,
subplot(221);semilogy(s),title('(a) singular values');
subplot(222);semilogy(f(nt/2+1:nt),fftshift(WW(nt/2+1:nt)));
title('(b) Spectrum of the wavelet');ylabel(ylab);xlabel('f(HZ)');
subplot(223);plot(U(:,3:97:100));title('(c) singular vector u_3 and u_{100}');
subplot(224);plot(V(:,3:97:100));title('(d) singular vector v_3 and v_{100}');
figure
subplot(311);plot(U(:,3:97:100));title('(a) singular vector u_3 and u_{100}');
subplot(312); picard(U,s,b);title('(b) clean data');
subplot(313); picard(U,s,d);title('(c) noisy data');

% Regularization parameters
% Use the L-curve criterion and GCV to determine the regularization
% parameters for Tikhonov regularization and truncated SVD.

figure
subplot(221);k_l1 = l_curve(U,s,b,'tsvd'); title('(a)'); %axis([1e-3,1,1,1e3]),      pause
subplot(222);k_gcv1 = gcv(U,s,b,'tsvd'); title('(b)');   %axis([0,20,1e-9,1e-1]),    pause
subplot(223);k_l2 = l_curve(U,s,d,'tsvd'); title('(c)'); %axis([1e-3,1,1,1e3]),      pause
subplot(224);k_gcv2 = gcv(U,s,d,'tsvd');  title('(d)');  %axis([0,20,1e-9,1e-1]),    pause

x_tsvd_gcv1 = tsvd(U,s,V,b,k_gcv1);
mytext1=sprintf('(b) Solution for exact data: SVD and %d singular vectors',k_gcv1);
x_tsvd_gcv2 = tsvd(U,s,V,d,k_gcv2);
mytext2=sprintf('(c) Solution for noisy data: SVD and %d singular vectors',k_gcv2);

figure,
subplot(311);plot(x);title('(a) Initial model');
subplot(312);plot(x_tsvd_gcv1);title(mytext1);
subplot(313);plot(x_tsvd_gcv2);title(mytext2);


% Iterative solution with prior information
% choose exact data or data with noise
b=d;  

nfig1=10
nfig2=11

x_k=x_tsvd_gcv2;
L=diag(1./abs(x_k));

niter=4;
figure
I = 1;
 
for i=1:4
  isvd=2*i;
  subplot(2,2,i); plot(1:nt,V(:,isvd)); axis([1,nt,-1,1])
  xlabel(['i = ',num2str(isvd)]), I = I + 1;
end
subplot(221);title('(a)')
subplot(222);title('(b)')
subplot(223);title('(c)')
subplot(224);title('(d)')

% TSVD solution
figure(nfig1),
subplot(211); k_tsvd = gcv(U,s,b,'tsvd');title('(a)');
X_I = tsvd(U,s,V,b,1:k_tsvd);
figure(nfig2)
subplot(3,1,1);
plot(1:n,X_I,1:n,x,'x'), axis([1,n,-1.2,1.2]), xlabel('L = I')
title('Zero order regularization');
axis1=axis;    

ifig2=1;
  

for iter=1:niter
  L=diag(1./(max(abs(x_k),eps1)));
  [UU,sm,XX] = cgsvd(A,L);
  axisvector=[1,nt,min(XX(:,nt)),max(XX(:,nt))];
  if (iter==niter)
    % Plot Generalized Eigenvectors
    figure
    I=1;  
    for i=[nt,nt-1,nt-2,nt-3]      
      subplot(2,2,I); plot(1:nt,XX(:,i)), 
      axis(axisvector);
      xlabel(['i = ',num2str(i)]), I = I + 1;
    end
    subplot(221);title('(a)')
    subplot(222);title('(b)')
    subplot(223);title('(c)')
    subplot(224);title('(d)')   
    %text(10,1.2,'GSVD: Right generalized singular vectors XX(:,i)')
  end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot GCV 

  figure(nfig1),
  subplot(2,1,2); k_tgsvd = gcv(UU,sm,b,'tsvd');  title('(b)');

  
  X_L = tgsvd(UU,sm,XX,b,1:k_tgsvd);
  % Plot GCV 
  if (iter>2)  
    figure(nfig2)
    ifig2=ifig2+1;
    subplot(3,1,ifig2);
    mytext=sprintf('iter=%d\n',iter);
    plot(1:n,X_L,1:n,x,'x'), axis([1,n,-1.2,1.2]),ylabel(mytext); 
    axis(axis1);
  end  
  x_k=X_L(:,k_tgsvd); 

end;
xsol=x_k;
figure(nfig2)
subplot(311);title('(a)')
subplot(312);title('(b)')
subplot(313);title('(c)')

figure,
subplot(311),semilogy(1:n,sm(:,1)./sm(:,2));title('(a) Generalized eigenvalues \gamma_i')
subplot(312),semilogy(1:n,sm(:,1));title('(b) \sigma_i')
subplot(313),semilogy(1:n,sm(:,2));title('(c) Generalized \mu_i')

figure,
subplot(211);plot(x);title('(a) Initial model');
subplot(212);plot(xsol);
mytext=sprintf('(b) Final solution with GSVD and %d eigenvectors\n',k_tgsvd);
title(mytext);









