%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOSC 555 - Problem set #3.                   
%
% Minimum of Rosenbrock function by using Newton's method 
% The search is performed by the function newton2dim
% The function, first and secon derivatives are in fx, d1fx, d2fx,
% etc. Also derivatives in parametric form are supplied to apply
% the line search
% 
% Exercise 1 is Newton method with line search
% Exercise 2 is Newton method without line search
% The minimum is check with norm(gradient) < tol
% and a evaluating the function at random points around the minimum
% (function check minimum)
%
% Daniel Trad    January 2001- EOSC 555  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all
exercise=1;
if (exercise==1) ex='Exercise 1'
else ex='Exercise 2'
end  

% Functions
fx='rosenbrock'  % f(x) function
d1fx=['d1'fx];   %  first derivative (gradient)
d2fx=['d2'fx];   %  second derivative (Hessian)
fa='rosenbrock_alpha'; % function in parametric form
d1fa=['d1'fa]; % First derivative of the function in parametric form
d2fa=['d2'fa]; % Second derivative of the function in parametric form
fx2='rosenbrock2'  % Function for steep.m that return f(x) and f'(x)

xg=-2:.2:2;
yg=-4:.2:3;
[xx,yy]=meshgrid(xg,yg);
zz=100*(yy-xx.^2).^2+(1-xx).^2;

figure(1); mesh(xg,yg,zz); xlabel('x');ylabel('y');title(['Our old friend:' ...
		    ' Rosenbrock']);

x0=[-1.2 1]; maxiter=230; tol=1e-3;

if (exercise==1)
  search=1; % Use line search
else
  search=0; % Do not use line search
end

[y,fy,normdeltay,normgrad,eigen,theta]=newton2dim(fx,d1fx,d2fx,fa,d1fa,d2fa,x0,maxiter,tol,search); % Newton

if (checkminimum(fx,[y(end,1) y(end,2)],0.5,tol)==0) 
  check='Minimum check OK'
else
  check='This is not a local minimum'
end

iter=1:max(size(y));
figure(2);contour(xg,yg,zz,20);hold on

i=1:iter(end);
plot(y(i,1),y(i,2));plot(y(end,1),y(end,2),'g*');
plot(y(i,1),y(i,2),'.');
hold off

text(x0(1),x0(2),'<--I started here!!!');
if (exercise==1) text(y(end,1),y(end,2),'<--with line search!!!');
else text(y(end,1),y(end,2),'<--WITHOUT line search, Jah!!!');
end

text(0,-1,check); 

text2=sprintf('xmin=[%4.2f %4.2f] Found in %d iterations \n',y(end,1),y(end,2), ...
	      iter(end));
text(-1.,1.5,text2);

if (exercise==1)    title('Minimization using Newton with line search')
else  title('Minimization using Newton without line search')
end
  
figure(3);
subplot(311);semilogy(fy),title(['f(x):',fx,'-',ex])
subplot(312);semilogy(normdeltay),title('\delta x')
subplot(313);semilogy(normgrad),title('|\Delta f|'); xlabel('iterations');

figure(4);
if (exercise==1)
  subplot(311);semilogy(eigen(:,1));
  title('Maximum eigenvalues for Mr. Hessian: With line search')
  subplot(312);semilogy(eigen(:,2));
  title('Minimum eigenvalues for Mr. Hessian: With line search')
else
  subplot(311);semilogy(eigen(:,1));
  title('Maximum eigenvalues for Mr. Hessian: No line search')
  subplot(312);semilogy(eigen(:,2));
  title('Minimum eigenvalues for Mr. Hessian: No line search')
end

I=find(eigen<=0.1)

if (length(I)==0) text(2,max(eigen(:,2))/2,'Mr Hessian has always been a SPD guy');
else text(2,max(eigen(:,2))/2,'Mr Hessian has not always been a SPD guy');
end


subplot(313); plot(theta),axis([0 length(theta)+1 40 90])
title('\theta between the Newton direction and the steepest descent  direction') 
xlabel('iterations');ylabel('\theta');



