% EOSC 555 --- Daniel Trad 
% Problem Set 1
close all
clear
x0=2.0;
x=[0 2]
tol=1e-3;
fx0='expx2'
fx='expx2minus3';
d1fx=['d1',fx0];
d2fx=['d2',fx0];
maxiter=40
z=0
v=[x(1) x(2) -3 10]

figure(1)

subplot(411); fplot(fx,x);  
axis(v);
ylabel(fx),xlabel('x')

hold on; 
plot(x0,1,'b*');  
%plot([x(1) x(1)],[0 1],'r')
%plot([x(2) x(2)],[0 1],'r')

[y,z] = newton(fx,d1fx,x0,tol,maxiter); 

if (feval(fx,y)<=tol)
    text=sprintf('-->found a zero at x=%5.2f with f(x)=%5.2f\n',y,feval(fx,y));
else
  text=sprintf('It did not find a root\n');
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

v2=[1 length(z) -10 100]

itx=1:length(z);
subplot(412);
%text2='function'
plot(itx,fz,itx,fz,'*'),ylabel(fx),xlabel('iterations')
axis(v2);


subplot(413);
plot(itx,d1fz,itx,d1fz,'*'),ylabel(['d/dx ',fx]),xlabel('iterations');
axis(v2);

subplot(414);
plot(itx,d2fz,itx,d2fz,'*'),ylabel(['d/dx ',fx]),xlabel('iterations')
axis(v2);

%print -dps2 ps1fig2.ps



