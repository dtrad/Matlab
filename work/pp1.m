......Continuation from previous page 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same as before with conjugate directions preconditioned
% Find conjugate directions
[x,rho,niter,J,S]=conj_direc(A,b,tol,length(xmin),x0);
S=S(:,1:2);
lambda=svd(A);
p1=S(:,1);
p2=S(:,2);
p10=p1;
p20=p2;
% Orthonormalize the p vectors
p1=p1/sqrt(p1'*A*p1); 
p2=p2/sqrt(p2'*A*p2);
S=[p1 p2];
%Test for conjugacy and orthgonality
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1
mindescdir(A,b,S,x0,z,x1,x2,6,'D= Identity');
% Try to converge in only one iteration
D=S'*A*S; bb=S'*b;
x0t=inv(S)*x0;
g=(D*x0t-bb);
figure(6)
subplot(224)
zz=1/2*(x1.*(D(1,1).*x1+D(1,2).*x2)+x2.*(D(2,1).*x1+D(2,2).*x2))-bb(1).*x1-bb(2).*x2;
contour(x1,x2,zz,20); 
hold on
ysd=desc_dir(D,bb,-g,x0t,maxiter,tol)
contour(x1,x2,zz,20); hold on
i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');
title('Steepest descent in rotated system')
hold off
return;
