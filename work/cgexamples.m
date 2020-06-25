% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000


A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1];
m=[1;2;3;4;5;0]
b=A*m;
[nn mm]=size(A);
m0=zeros(mm,1);
W=diag([1 1 1 1 1]);   % Data Preconditioner 
M=diag([1 1 1 1 1 1]); % Model preconditioner
[U,S,V]=svd(A);
v6=V(:,6); % Null space
tol=1e-10;

% Perfect data, underdetermined problem, zero null space.

[x,rho,niter]=lpcgnr(A,b,M,tol,m0);
[m x]
rho

[x,rho,niter]=mpcgnr(A,b,W,tol,m0);
[m x]
rho

[x,rho,niter]=wpcgnr(A,b,M,tol,m0,5,W);
[m x]
rho

[x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
[m x]
rho

%pause;
 
% Use prior information
M=diag([1 1 1 1 1 1e5]);

[x,rho,niter]=lpcgnr(A,b,M,tol,m0);
[m x]
rho

[x,rho,niter]=mpcgnr(A,b,W,tol,m0);
[m x]
rho

[x,rho,niter]=wpcgnr(A,b,M,W,tol,m0);
[m x]
rho

[x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
[m x]
rho

%pause

% Data error, use data weight
%display('Data error, use data weight and prior information');
b0=b;
b(5)=10;
M=diag([1 1 1 1 1 1e5]);
W=diag([1 1 1 1 2]);

[x,rho,niter]=lpcgnr(A,b,M,tol,m0);
[m x] 
[b A*x]
rho

[x,rho,niter]=mpcgnr(A,b,W,tol,m0);
[m x ]
rho
[b0 b A*x]

[x,rho,niter]=wpcgnr(A,b,M,W,tol,m0);
[m x]
rho
[b0 b A*x]

[x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
[m x]
rho
[b0 b A*x]
%pause
% starting model with Null space
b=b0;
W=diag([1 1 1 1 1]);
M=diag([1 1 1 1 1 1]);

[x,rho,niter]=lpcgnr(A,b,M,tol,m0);
for i=1:6, rhox(i)=x'*V(:,i),end
[(m-x)./v6]
rho

[x,rho,niter]=mpcgnr(A,b,W,tol,m0);
[(m-x)./v6]
rho

[x,rho,niter]=wpcgnr(A,b,M,W,tol,m0);
[(m-x)./v6]
rho

[x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
[(m-x)./v6]
rho
%pause

%display('Now with prior information');
W=diag([1 1 1 1 1]);
M=diag([1 1 1 1 1 1e5]);

[x,rho,niter]=lpcgnr(A,b,M,tol,m0);
[(m-x)./v6]
rho

[x,rho,niter]=mpcgnr(A,b,W,tol,m0);
[(m-x)./v6]
rho

[x,rho,niter]=wpcgnr(A,b,M,W,tol,m0);
[(m-x)./v6]
rho

[x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
[(m-x)./v6]
rho
%pause