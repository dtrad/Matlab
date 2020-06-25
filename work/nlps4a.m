%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOSC 555 - Problem set #4.                   
% Exercise 1
% CG solver for solving a system of equations for PD matrix
% Exercise 2
% Minimum of a quadratic function by using steepest descent and
% full Newton's method 
% The search is performed by the function steep_desc2 using exact
% line search.
% The function, first and second derivatives are in fx, d1fx,H
% 
% The minimum is check with norm(gradient) < tol
%
% Daniel Trad    February 2001- EOSC 555  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing Conjugate gradients algorithm
% Daniel Trad - UBC March-15 2000

clear
close all
exercise=2;
  
%exercise 1
if (exercise==1)
  A=[3 5 5 ; 5 13 9 ; 5 9 11];
  b=[13 23 29];b=b(:);
  tol=1e-10;
  [x,rho,niter]=cg555(A,b,tol);
  'b A*x'
  [b A*x]
  [rho' x ] 
  niter
else
  fx='funcps4';
  d1fx=['d1',fx];
  xg=-4:.2:4;
  yg=-4:.2:4;
  [xx,yy]=meshgrid(xg,yg);
  zz=4*xx.^2+3*yy.^2-4*xx.*yy+xx+10;
  
  % Starting point
  x0=[2;3];
  % Hessian  
  H=[8 -4; -4 6];
  tol=1e-5;  
  maxiter=100;
  % Steepest descent  
  [ysd,fy,normdeltay,normgrad]=steep_desc2(fx,d1fx,H,x0,maxiter,tol); 
  figure(1);
  contour(xg,yg,zz,20); hold on
  i=1:max(size(ysd));
  plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
  plot(ysd(i,1),ysd(i,2),'r.');
  title('Steepest descent for quadratic functions and full Newton');
  text1=sprintf('Steepest descent need %d iterations',i(end));
  text(0.5,0,text1);
  xlabel('x');ylabel('y');
  figure(2);
  subplot(311);plot(fy,'o'),title('f(x):')
  subplot(312);semilogy(normdeltay,'o'),title('||\delta x||')
  subplot(313);semilogy(normgrad,'o'),title('||\nabla f||'); xlabel('iterations');
  
  figure(1)
  % Full Newton step
  % Gradient
  g=feval(d1fx,x0);
  [x,rho,niter]=cg555(H,-g,tol);
  '[-g H*x]'
  [-g H*x]
  '[rho x]'
  [rho' x]
  display(sprintf('niter=%d',niter))
  y=x0+x;
  
  plot(x0(1),x0(2),'b+');plot(y(1),y(2),'g*');
  plot([x0(1) y(1)],[x0(2) y(2)]);
  hold off
  theta=(180/pi)*acos(dot(x,-g)/(norm(x)*norm(-g)));
  text2=sprintf('Angle Newton-initial steep. dec. %5.2f degrees ',theta);
  text(-3,2,text2);
end




