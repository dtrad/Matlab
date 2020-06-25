% Testing Conjugate gradients algorithms with model and data
% preconditioner. Influence on outliers, prior information, and
% null space.
% Daniel Trad - UBC March-15 2000
echo on

A=[1 2 3 4 5 1; 1 0 -1 0 2 1;2 5 2 1 1 1;1 1 1 1 1 1;1 -1 -1 1 1 1;2 ...
  3 1 3 4 1; 1 1 1 7 2 1];
m=[1;2;3;4;5;0]
b=A*m;
[nn mm]=size(A);
m0=zeros(mm,1);
W=diag([1 1 1 1 1 1 1]);   % Data Preconditioner 
M=diag([1 1 1 1 1 1]); % Model preconditioner
[U,S,V]=svd(A);
v6=V(:,6); % Null space
tol=1e-10;

test=2


switch test
  % Perfect data, underdetermined problem, zero null space.
 case 1,
  [x,rho,niter]=lpcgnr(A,b,M,tol,m0);
  [m x]
  rho
  
  [x,rho,niter]=mpcgnr(A,b,W,tol,m0);
  [m x]
  rho

  display('test wpcgnr');
  
  [x,rho,niter]=wpcgnr(A,b,M,tol,m0,length(m0),W);
  [m x]
  rho

  [x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
  [m x]
  rho

%pause;
 case 2, 
  % Use prior information
  b0=b;
  %b(2)=b(2)+1;
  M=diag([1 1 1 1 1 .1]);
  % Make the problem underdetermined by zeroing the weights for 
  % some of the equations
  % The problem becomes more underdetermined and the prior
  % information in the preconditioner M is more important
  % The perfect preconditioner is when m_true is 
  W=diag([1 0 0 1 0 1 1]); 
  M=diag(m+0.01);
  itercg=6;
  if (0)
    [x,rho,niter]=lpcgnr(A,b,M,tol,m0);
    [m x] 
    [b A*x]
    rho
  elseif(0)
    [x,rho,niter]=mpcgnr(A,b,W,tol,m0);
    [m x ]
    rho
    [b0 b A*x]
  elseif(1) 
    [x,rho,niter]=wpcgnr(A,b,M,tol,m0,itercg,W);
    [x2,rho2,niter2]=wpcgnr(A,b,inv(M),tol,m0,itercg,W);
    [m x x2]
    [rho(:) rho2(:)]
    [b0 b A*x A*x2]
  elseif(0)
    [x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
    [m x]
    rho
    [b0 b A*x]
  end
  %pause

 case 3,
  % Data error, use data weight
  %display('Data error, use data weight and prior information');
  b0=b;
  b(5)=100;
  M=diag([1 1 1 1 1 1]);
  W=diag([1 1 1 1 10 1 1]);
  itercg=6;
  if (0)
  [x,rho,niter]=lpcgnr(A,b,M,tol,m0);
  [m x] 
  [b A*x]
  rho
  elseif(0)
  [x,rho,niter]=mpcgnr(A,b,W,tol,m0);
  [m x ]
  rho
  [b0 b A*x]
  elseif(1) 
    [x,rho,niter]=wpcgnr(A,b,M,tol,m0,itercg,W);
    [x2,rho2,niter2]=wpcgnr(A,b,M,tol,m0,itercg,inv(W));
    % the last one is the right. Use small values to penalize outliers)
    [m x x2]
    [rho(:) rho2(:)]
    [b0 b A*x A*x2]
  elseif(0)
  [x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
  [m x]
  rho
  [b0 b A*x]
  end
  %pause
 case 4,
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
  
  [x,rho,niter]=wpcgnr(A,b,M,tol,m0,length(m0),W);
  [(m-x)./v6]
  rho
  
  [x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
  [(m-x)./v6]
  rho
  %pause
 case 5,
  %display('Now with prior information');
  W=diag([1 1 1 1 1]);
  M=diag([1 1 1 1 1 1e5]);
  
  [x,rho,niter]=lpcgnr(A,b,M,tol,m0);
  [(m-x)./v6]
  rho
  
  [x,rho,niter]=mpcgnr(A,b,W,tol,m0);
  [(m-x)./v6]
  rho
  
  [x,rho,niter]=wpcgnr(A,b,M,tol,m0,length(m0),W);
  [(m-x)./v6]
  rho
  
  [x,rho,eta,niter]=wtcgls(A,M.^2,b,5,tol,1); 
  [(m-x)./v6]
  rho
  %pause
 otherwise, 
  display('no test performed............');
end
