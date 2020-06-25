%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOSC 555 - Problem set #5.    Exercise 1
% Minimum of a quadratic function by using 
% a) steepest descent,
% b) Eigenvector directions
% c) CG directions
% d) Conjugate directions from Gram Schmidt 
% Daniel Trad    February 2001- EOSC 555  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
clear
A=[3 2;2 6]
H=A; % Hessian
b=[2;-8];
fx='funcps5';
d1fx=['d1',fx];
tol=1e-5;  
maxiter=40;
[U,S,V]=svd(A);
xg=-5:.2:5;
yg=-5:.2:5;
[x1,x2]=meshgrid(xg,yg);
z=1/2*(x1.*(A(1,1).*x1+A(1,2).*x2)+x2.*(A(2,1).*x1+A(2,2).*x2))-b(1).*x1-b(2).*x2;
contour(x1,x2,z,20); hold on
xmin=A\b;
plot(xmin(1),xmin(2),'g*');
x0=[-2;2];
plotvector(V(:,1),xmin,'g');
plotvector(V(:,2),xmin,'b');
text(2.2,-1.5,'<---Eigenvector 1')
text(2.5,-2.3,'<---Eigenvector 2')
xaxis=[1 0]
theta1=(180/pi)*acos(dot(xaxis,V(:,1))/(norm(xaxis)*norm(V(:,1))))
theta2=(180/pi)*acos(dot(xaxis,V(:,2))/(norm(xaxis)*norm(V(:,2))))

% Steepest Descent minimization
[ysd,fy,normdeltay,normgrad]=steep_desc2(fx,d1fx,H,x0,maxiter,tol); 
i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');
text(-4,0,'Steepest Descent--->');
title('Quadratic Function');
xlabel('x');
ylabel('y');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimization along eigenvectors
% Test for eigenvectors
p1=V(:,1);
p2=V(:,2);
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1

S=[p1 p2];
mindescdir(A,b,S,x0,z,x1,x2,2,'Eigenvectors');
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same as before with conjugate directions
% Find conjugate directions
[x,rho,niter,J,S]=conj_direc(A,b,tol,length(x0),x0);
S=S(:,1:2);
p1=S(:,1);
p2=S(:,2);

%Test for conjugacy
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1
mindescdir(A,b,S,x0,z,x1,x2,3,'Conjugate Gradient');

display('The same as before with Grand Schmidt') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find conjugate directions using Gram-Schmidt
% [x,rho,niter,J,S]=conj_direc(A,b,tol,length(xmin),[0 0]);
% Case 1: Firs vecotr is minus gradient
theta=pi/2;
rot=[cos(theta) -sin(theta); sin(theta) cos(theta)];
u1=A*x0-b;
u2=rot*u1; % Second vetor is perpendicular to first one
% Grand Schmidt
p1=u1/norm(u1);
p2=u2-(u2'*A*p1)/(p1'*A*p1)*p1
p2=p2/norm(p2);
S=[p1 p2];
display('Test for Conjugacy') 
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1
mindescdir(A,b,S,x0,z,x1,x2,4,['Grand Schmidt: case 1, first direction' ...
	   'is the gradient ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 2: Firs vector is x axis
u1=[1;0]
u2=rot*u1; % Second vetor is perpendicular to first one
% Grand Schmidt
p1=u1/norm(u1);
p2=u2-(u2'*A*p1)/(p1'*A*p1)*p1
p2=p2/norm(p2);
S=[p1 p2];
display('Test for Conjugacy') 
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1
mindescdir(A,b,S,x0,z,x1,x2,5,['Gram Schmidt: case 2 first direction' ...
	   'is the x axis']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

it continues on next page



