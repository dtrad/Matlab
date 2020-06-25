function [xsol]=regudecon(A,b_bar,x,L,e,neig,niter,eps1,eps2)
% Christian Hansen, IMM, 12/19/97.
% Adapted for deconvolution.
% Daniel Trad- UBC - 13-07-99 
echo on, 

% Part 1.  The discrete Picard condition
% --------------------------------------
%
% First generate a "pure" test problem where only rounding
% errors are present.  Then generate another "noisy" test
% problem by adding white noise to the right-hand side.
%
% Next compute the SVD of the coefficient matrix A.
%
% Finally, check the Picard condition for both test problems
% graphically.  Notice that for both problems the condition is
% indeed satisfied for the coefficients corresponding to the
% larger singular values, while the noise eventually starts to
% dominate.
% 


if (nargin<5)
     randn('seed',41997);
     e = 1e-2*randn(size(b_bar)); 
end

if (nargin<6) neig=5;end
if (nargin<7) niter=4;end
if (nargin<8) eps1=1e-10;end
if (nargin<9) eps2=1e-10;end
b = b_bar + e;
figure,
subplot(211);plot(b);title('data with noise');
subplot(212);plot(b_bar);title('data without noise');
%b_bar=zeros(size(b));

[U,s,V] = csvd(A);
figure,
subplot(2,1,1); picard(U,s,b_bar);
subplot(2,1,2); picard(U,s,b);


% Part 2.  Filter factors
% -----------------------
%
% Compute regularized solutions to the "noisy" problem from Part 1 
% by means of Tikhonov s method and LSQR without reorthogonalization.
% Also, compute the corresponding filter factors.
%
% A surface (or mesh) plot of the solutions clearly shows their dependence
% on the regularization parameter (lambda or the iteration number).

lambda = [3e1,1e1,3,1,3e-1,1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,3e-6,1e-6];
X_tikh = tikhonov(U,s,V,b,lambda);
F_tikh = fil_fac(s,lambda);
iter = 30; reorth = 0;
[X_lsqr,rho,eta,F_lsqr] = lsqr(A,b,iter,reorth,s);
figure
subplot(2,2,1); surf(X_tikh), title('Tikhonov solutions'), axis('ij')
subplot(2,2,2); surf(log10(F_tikh)), axis('ij')
                title('Tikh filter factors, log scale')
subplot(2,2,3); surf(X_lsqr(:,1:17)), title('LSQR solutions'), axis('ij')
%subplot(2,2,4); surf(log10(F_lsqr(:,1:17))), axis('ij')
title('LSQR filter factors, log scale')


% Part 3.  The L-curve
% --------------------
%
% Plot the L-curves for Tikhonov regularization and for
% LSQR for the "noisy" test problem from Part 1.
%
% Notice the similarity between the two L-curves and thus,
% in turn, by the two methods.
figure
subplot(1,2,1); l_curve(U,s,b); %axis([1e-3,1,1,1e3])
subplot(1,2,2); plot_lc(rho,eta,'o'); %axis([1e-3,1,1,1e3])


% Part 4.  Regularization parameters
% ----------------------------------
%
% Use the L-curve criterion and GCV to determine the regularization
% parameters for Tikhonov regularization and truncated SVD.
%
% Then compute the relative errors for the four solutions.
figure
subplot(221);lambda_l = l_curve(U,s,b);   %axis([1e-3,1,1,1e3]),      pause
subplot(222);k_l = l_curve(U,s,b,'tsvd'); %axis([1e-3,1,1,1e3]),      pause
subplot(223);lambda_gcv = gcv(U,s,b);     %axis([1e-6,1,1e-9,1e-1]),  pause
subplot(224);k_gcv = gcv(U,s,b,'tsvd');   %axis([0,20,1e-9,1e-1]),    pause

%lambda_l=2;
k_l=k_gcv;

x_tikh_l   = tikhonov(U,s,V,b,lambda_l);
x_tikh_gcv = tikhonov(U,s,V,b,lambda_gcv);
if isnan(k_l)
  x_tsvd_l = zeros(size(x)); % Spline Toolbox not available.
else
  x_tsvd_l = tsvd(U,s,V,b,k_l);
end
x_tsvd_gcv = tsvd(U,s,V,b,k_gcv);
results=[norm(x-x_tikh_l),norm(x-x_tikh_gcv),...
 norm(x-x_tsvd_l),norm(x-x_tsvd_gcv)]/norm(x)
%II=find(results==min(results-1))
II=input('According to these results which method for prior inf: 1,2,3: ');

figure,
subplot(221);plot(x_tikh_l);
title(sprintf('Tikhonov with lambda=%f\n',lambda_l));
subplot(222);plot(x_tikh_gcv);
title(sprintf('Tikhonov with GVC k=%d\n',lambda_gcv));
subplot(223);plot(x_tsvd_l);
title(sprintf('Final solution with TSVD and %d eigenvectors\n',k_l));
subplot(224);plot(x_tsvd_gcv);
title(sprintf('tsvd with GCV k=%d \n',k_l));

if (II==1) xp=x_tikh_l;
elseif (II==2) xp=x_tikh_gcv;
elseif (II==3) xp=x_tsvd_gcv;
end;

%xp=x_tsvd_gcv;
%xp=x_tikh_gcv; 

% Part 5.  Standard form versus general form
% ------------------------------------------
%
% Generate a new test problem: inverse Laplace transformation
% with white noise in the right-hand side.
%
% For the general-form regularization, choose minimization of
% the first derivative.
%
% First display some left singular vectors of SVD and GSVD; then
% compare truncated SVD solutions with truncated GSVD solutions.
% Notice that TSVD cannot reproduce the asymptotic part of the
% solution in the right part of the figure.

n = length(x); 
%[A,b,x] = ilaplace(n,2);
%b = b + 1e-4*randn(size(b));
%L = get_l(n,1);
niter=4;

for iter=1:niter
 if (iter>1) xp=xsol;end
 L=diag(1./(abs(xp)+eps2)+eps1);
 [UU,sm,XX] = cgsvd(A,L);
 if (iter==1) [U,s,V] = csvd(A); 

 figure,
 subplot(311),semilogy(1:n,sm(:,1)./sm(:,2));title('Generalized eigenvalues \gamma_i')
 subplot(312),semilogy(1:n,sm(:,1));title('\sigma_i')
 subplot(313),semilogy(1:n,sm(:,2));title('Generalized \mu_i')
 %subplot(213),semilogy(1:n,sm(:,1).^2+sm(:,2).^2));


 figure
 I = 1;
 for i=[1,2,3,4]
  subplot(2,2,I); plot(1:n,V(:,i)); axis([1,n,-1,1])
  xlabel(['i = ',num2str(i)]), I = I + 1;
 end
 subplot(2,2,1), text(12,1.2,'SVD: Right singular vectors V(:,i)'), 
end;


figure
I = 1;
for i=[n,n-1,n-2,n-3]
  subplot(2,2,I); plot(1:n,XX(:,i)), axis([1,n,-1,1]);
  xlabel(['i = ',num2str(i)]), I = I + 1;
end
subplot(2,2,1)
text(10,1.2,'GSVD: Right generalized singular vectors XX(:,i)')
figure,
subplot(211); k_tsvd = gcv(U,s,b,'tsvd')
subplot(212); k_tgsvd = gcv(UU,sm,b,'tsvd') 


%k_tsvd = 7; k_tgsvd = 6;
if (iter==1) X_I = tsvd(U,s,V,b,1:k_tsvd);end;
X_L = tgsvd(UU,sm,XX,b,1:k_tgsvd);
figure
subplot(2,1,1);
  plot(1:n,X_I,1:n,x,'x'), axis([1,n,-1.2,1.2]), xlabel('L = I')
axis1=axis;
subplot(2,1,2);
plot(1:n,X_L,1:n,x,'x'), axis([1,n,-1.2,1.2]), xlabel('L \neq I')
axis(axis1);
xsol=X_L(:,k_tgsvd); 

figure,
subplot(211);plot(x);title('Initial model');
subplot(212);plot(xsol);
mytext=sprintf('Final solution with GSVD and %d eigenvectors\n',k_tgsvd);
title(mytext);

end

% Part 6.  No square integrable solution
% --------------------------------------
%
% In the last example there is no square integrable solution to
% the underlying integral equation (NB: no noise is added).
%
% Notice that the discrete Picard condition does not seem to
% be satisfied, which indicates trouble!

%[A,b] = ursell(32); [U,s,V] = csvd(A);
%picard(U,s,b); pause

% This concludes the demo.
echo off











