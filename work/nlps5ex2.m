%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOSC 555 - Problem set #5.                   
%
% Minimum of Rosenbrock function by using nonlinear CG method 
% The CG search is in the function nlcgps5.m
% The linear search is performed by the function newtonparam
% The function, first and secon derivatives are in fx, d1fx, d2fx,
% are given.
% The minimum is check with norm(gradient) < tol
%
% Daniel Trad    Febrary 2001- EOSC 555  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all
ex='Exercise 2'

% Functions
fx='rosenbrock'  % f(x) function
d1fx=['d1'fx];   %  first derivative (gradient)
d2fx=['d2'fx];   %  second derivative (Hessian)

xg=-2:.2:2;
yg=-4:.2:3;
[xx,yy]=meshgrid(xg,yg);
zz=100*(yy-xx.^2).^2+(1-xx).^2;

figure(1); mesh(xg,yg,zz); xlabel('x');ylabel('y');title(['Our old friend:' ...
		    ' Rosenbrock']);

x0=[-1.2 1]; maxiter=1000; tol=1e-3;

[y,fy,normdeltay,normgrad]=nlcgps5(fx,d1fx,d2fx,x0,maxiter,tol);

iter=1:max(size(y));
figure(2);contour(xg,yg,zz,100);hold on

i=1:iter(end);
plot(y(i,1),y(i,2));plot(y(end,1),y(end,2),'g*');
plot(y(i,1),y(i,2),'.');
hold off


text(x0(1),x0(2),'<--I started here!!!');
text(y(end,1),y(end,2),'<--Minimum!!!');
text2=sprintf('xmin=[%4.2f %4.2f] Found in %d iterations \n',y(end,1),y(end,2), ...
	      iter(end));
text(-1.,1.5,text2);

title('Minimization using Noninear CG')

  
figure(3);
subplot(311);semilogy(fy),title(['f(x):',fx,'-',ex])
subplot(312);semilogy(normdeltay),title('\Delta x')
subplot(313);semilogy(normgrad),title('|\nabla f|'); xlabel('iterations');

return;

