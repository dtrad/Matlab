%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOSC 555 - Problem set #2.                   Daniel Trad    %
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all
exercise=2;
if (exercise==1) ex='Exercise 1'
else ex='Exercise 2'
end  

% Functions
fx='rosenbrock'  % f(x) function
d1fx=['d1'fx];   %  first derivative
fa='rosenbrock_alpha'; % function in parametric form
d1fa=['d1'fa]; % First derivative of the function in parametric form
d2fa=['d2'fa]; % Second derivative of the function in parametric form
fx2='rosenbrock2'  % Function for steep.m that return f(x) and f'(x)

xg=-2:.2:2;
yg=-1:.2:3;
[xx,yy]=meshgrid(xg,yg);
zz=100*(yy-xx.^2).^2+(1-xx).^2;

figure(1); mesh(xg,yg,zz); xlabel('x');ylabel('y');title('Rosenbrock');

x0=[-1.2 1]; maxiter=230; tol=1e-3;

if (exercise==1)
  [y,fy,normdeltay,normgrad]=steep_desc(fx,d1fx,fa,d1fa,d2fa,x0, ...
					maxiter,tol); % My steep descent
else
  [x,histout,costdata,y]=steep(x0(:),fx2,tol,maxiter);  % Kelley
  normgrad=histout(:,1);
  fy=histout(:,2);
  deltay=y(2:end,:)-y(1:end-1,:);
  deltay=[deltay;zeros(size(y(1,:)))];
  for ii=1:max(size(deltay)) 
    normdeltay(ii)=norm(deltay(ii,:));
  end
end	  

iter=1:max(size(y));
figure(2);contour(xg,yg,zz,20);hold on

i=1:iter(end);
plot(y(i,1),y(i,2));plot(y(end,1),y(end,2),'g*');hold off

text(x0(1),x0(2),'<--I started here!!!');
if (exercise==1) text(y(end,1),y(end,2),'<--I found you!!!');
else text(y(end,1),y(end,2),'<--Mr Kelley found you!!!');
end
text2=sprintf('xmin=[%4.2f %4.2f] Found in %d iterations \n',y(end,1),y(end,2), ...
	      iter(end));
text(-1.,1.5,text2);

if (exercise==1)    title('Minimization using my steepest descent')
else  title('Minimization using steep from Kelley')
end
  
figure(3);first=230;
subplot(311);semilogy(fy(1:first)),title(['f(x):',fx,'-',ex])
subplot(312);semilogy(normdeltay(1:first)),title('\delta x')
subplot(313);semilogy(normgrad(1:first)),title('|\Delta f|'); xlabel('iterations');

