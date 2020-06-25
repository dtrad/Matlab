% Solve the deconvolution problem using svd for circulant matrices
% Daniel Trad - April 5 2000


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
freq=0.1;
dt=1;
nt=128;n=nt;
t=0:nt-1;t=t*dt;
niter=4;

% 
% Prior information (Cm^{-1/2}=1/(abs(x)+eps2)+eps1;

eps1=1e-5;  
eps2=1e-5;

w=rickerm(freq,dt);
nw=length(w);
%ww=[w((nw+1)/2:nw)];
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
ww2=reverse(ww);
ww3=[ww2(nt) ww2(1:nt-1)];

% kernel
%A=toeplitz(ww,ww3);
%A=A(1:nt,1:nt);
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
subplot(221);plot(t(1:nw),w);title('Ricker wavelet');xlabel(xlab);ylabel(ylab)
subplot(222);plot(t,x);title('Reflectivity');xlabel(xlab);ylabel(ylab)
subplot(223);plot(t,b);title('Reflectivity+noise');xlabel(xlab);ylabel(ylab)
subplot(224);plot(t,d);title('solution');xlabel(xlab);ylabel(ylab)


% Let us studt first the kernel A and the data
% Singular vectors do not precise to be computed are just
% the fft coefficients of the first column.

ff=fft(A(:,1));
sf=(reverse(sort(abs((ff)))));
[U,ss,V]=svd(A);
figure,

semilogy(t,sf,'+',t,diag(ss),'o');
axis([0 t(nt) 1e-5 1e3])

% Now we create the singular vectors,

tt=0:nt-1;
j=0:nt-1;
W=exp(i*2*pi/nt*tt(:)*j(:).')/sqrt(nt);

for j=0:20:nt-1,
  figure
  %subplot(211),plot(real(exp(-i*tt*j*2*pi/nt))/nt),
  subplot(211),plot(real(W(:,j+1)));
  subplot(212),plot(V(:,j+1)),
end
figure
AA=W*diag(fft(A(:,1)))*W';

subplot(311),plot(real(A(:,1)));
xx=tsvd(W,ff,W,b,nt);
subplot(312),plot(b);
subplot(313),plot(tt,real(xx),tt,x,'+')

figure
for k=1:nt,
  xx=tsvd(W,ff,W,b,k);
  subplot(311),plot(t,real(xx),t,x,'+'),k,title('(a)');
  mytext=sprintf('%d eigenvetors',k);ylabel(mytext);
  figure(gcf),
end

for k=1:13,
  xx=tsvd(W,ff,W,d,k);
  subplot(312),plot(t,real(xx),t,x,'+'),k,title('(b)');
  mytext=sprintf('%d eigenvetors',k);ylabel(mytext);
  figure(gcf),
end
k=k+1;
xx=tsvd(W,ff,W,d,k);
subplot(313),plot(t,real(xx),t,x,'+'),k,title('(c)');
mytext=sprintf('%d eigenvetors',k);ylabel(mytext);



%LI=diag(abs(x)+1e-3);

%figure,
%for k=1:nt,xx=LI*tsvd(W,ff,W,d,k);
  %plot(t,real(xx),t,x,'+'),k,
  %figure(gcf),
%end