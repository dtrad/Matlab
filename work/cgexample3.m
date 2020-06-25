% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000

%Example 1
% Perfect data, underdetermined problem, zero null space.
example=3
if (example==1)
  A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1];
  m=[1;2;3;4;5;0]
  b=A*m;
  [nn mm]=size(A);

  W=diag([1 1 1 1 1]);   % Data Preconditioner 
  M=diag([1 1 1 1 1 1e-5]); % Model preconditioner

elseif (example==2) 
    
  A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1];
  
  m=[1;2;3;4;4;0]
  b=A*m;
  [nn mm]=size(A);
  m0=zeros(mm,1);
  W=diag([1 1 1 1 1]);   % Data Preconditioner 
  M=diag([1 1 1 1 1 1e-5]); % Model preconditioner
  tol=1e-15;

elseif (example==3)
  tol=1e-15
  small=1e-7;
  load input4.mat
  m=x;
  b=A*m;
  [nn mm]=size(A);
  W=diag(ones(nn,1));
  M=diag(ones(mm,1));
end
tol=1e-15;
m0=zeros(mm,1);
[x,rho,niter]=mpcgne(A,b,M,tol,m0);


[m x]
[b A*x]
format short e
rho
format 





