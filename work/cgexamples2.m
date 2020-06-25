% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000


A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1];
m=[1;2;3;4;5;0];
b=A*m;
[nn mm]=size(A);
m0=zeros(mm,1);
M=diag([1 1 1 1 1]);   % Data Preconditioner 
W=diag([1 1 1 1 1 1]); % Model preconditioner
[u,s,v]=svd(A'*A);
u0=u(:,6); % Null space
tol=1e-10;
iter=5;
b(2)=50;
[x,rho,niter]=cgne(A,b,M,W,tol,m0);
X=[m x];
for i=1:iter
  r=A*(m-x);
  M=diag([max(abs(r),1e-10)]);
  [x,rho,niter]=cgne(A,b,M,W,tol,m0);
  X=[X x];
end
X
