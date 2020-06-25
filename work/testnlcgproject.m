function [x]=testnlcgproject(lambda)
% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000

clean
close all


A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1];
m=[1;2;3;4;5;0];
b=A*m;
[nn mm]=size(A);
m0=zeros(mm,1);
M=diag([1 1 1 1 1]);   % Data Preconditioner 
W=diag([1 1 1 1 1 1e-3]); % Model preconditioner
[u,s,v]=svd(A'*A);
u0=u(:,6); % Null space
tol=1e-10;
iter=50;
%b(2)=50;
[xls,delnew,niter]=nlcgproject(A,b,tol,zeros(size(m)),W,0,iter)
close(1)
%lambda=1,
figure(1)
subplot(211);plot(m,'x');hold on
[xn,delnew,niter]=nlcgproject3(A,b,tol,zeros(size(m)),W,lambda,iter,0);
plot(xn,'r');axis([1 length(m) m(1) max(m)]);
hold off
subplot(212);plot(m,'x');
hold on
[xp,delnew,niter]=nlcgproject3(A,b,tol,zeros(size(m)),W,lambda,iter,1);
plot(xp,'r');axis([1 length(m) m(1) max(m)]);

hold off
X=[m xn xp  xls]
R=[b A*xn A*xp A*xls]

