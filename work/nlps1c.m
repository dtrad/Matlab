% EOSC 555 --- Daniel Trad 
% Problem Set 1
close all
clear
x0=3.5;
x=[0 6]
tol=1e-3;
fx='sin2';
d1fx=['d1',fx];
d2fx=['d2',fx];
maxiter=20
z=0

figure(1)

subplot(411); fplot(fx,x);  

ylabel(fx),xlabel('x')

hold on; 
plot(x0,1,'b*');  
%plot([x(1) x(1)],[0 1],'r')
%plot([x(2) x(2)],[0 1],'r')

[y,z] = newton(fx,d1fx,x0,tol,maxiter); 

if (feval(fx,y)<=tol)
    text=sprintf('-->found a zero at x=%5.2f with f(x)=%5.2f\n',y,feval(fx,y));
else
  text=sprintf('-->It did not find an zero\n');
end

title(['function: ',fx,'-starting point: ',num2str(x0),text]);
%while (feval(d1fx,y)<tol & feval(d2fx,y)<0)
%  x0=rand*x0;
%  [y,z] = newton(fx,d1fx,d2fx,x,x0,tol,maxiter); 
%end




iteraxis=(length(z)-1:-1:0)/length(z);
plot(z,iteraxis,'g*'); 
plot(y,0,'ro'); 
hold off;

for i=1:length(z)
  fz(i)=feval(fx,z(i));
end

for i=1:length(z)
  d1fz(i)=feval(d1fx,z(i));
end

for i=1:length(z)
  d2fz(i)=feval(d2fx,z(i));
end

itx=1:length(z);

subplot(412);
%text2='function'
plot(itx,fz,itx,fz,'*'),ylabel(fx),xlabel('iterations')



subplot(413);
plot(itx,d1fz,itx,d1fz,'*'),ylabel(['d/dx ',fx]),xlabel('iterations')


subplot(414);
plot(itx,d2fz,itx,d2fz,'*'),ylabel(['d/dx ',fx]),xlabel('iterations')


%print -dps2 ps1fig2.ps










