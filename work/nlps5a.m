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

% Steepest Descent minimization
[ysd,fy,normdeltay,normgrad]=steep_desc2(fx,d1fx,H,x0,maxiter,tol); 
i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off


% Minimization along eigenvectors
%Test for eigenvectors
test1=V(:,1)'*A*V(:,2)
test2=V(:,2)'*A*V(:,1)
test3=V(:,1)'*V(:,2)
test4=V(:,2)'*V(:,1)

p1=V(:,1);
p2=V(:,2);

S=[p1 p2];
SI=inv(S);

bb=S'*b;
D=S'*A*S;


%xx1=SI(1,1).*x1+SI(1,2).*x2;
%xx2=SI(2,1).*x1+SI(2,2).*x2;

figure(2)
subplot(221)
contour(x1,x2,z,20); hold on
[ysd,fy,normdeltay,normgrad]=desc_dir(A,b,S,[0 0],maxiter,tol)

i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off

subplot(222)
zz=1/2*(x1.*(D(1,1).*x1+D(1,2).*x2)+x2.*(D(2,1).*x1+D(2,2).*x2))-bb(1).*x1-bb(2).*x2;
contour(x1,x2,zz,20); hold on
ysdt=SI*ysd';
ysdt=ysdt';

i=1:max(size(ysd));
plot(ysdt(i,1),ysdt(i,2),'r');plot(ysdt(end,1),ysdt(end,2),'r*');
plot(ysdt(i,1),ysdt(i,2),'r.');

hold off


% The same as before with conjugate directions
% Find conjugate directions
[x,rho,niter,J,S]=conj_direc(A,b,tol,length(xmin),[0 0]);

S=S(:,1:2);

p1=S(:,1);
p2=S(:,2);


%Test for conjugacy
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1


bb=S'*b;
D=S'*A*S;
SI=inv(S);




subplot(223)
contour(x1,x2,z,20); hold on
[ysd,fy,normdeltay,normgrad]=desc_dir(A,b,S,[0 0],maxiter,tol)
i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off

subplot(224)
zz=1/2*(x1.*(D(1,1).*x1+D(1,2).*x2)+x2.*(D(2,1).*x1+D(2,2).*x2))-bb(1).*x1-bb(2).*x2;
contour(x1,x2,zz,20); hold on
ysdt=SI*ysd';
ysdt=ysdt';

i=1:max(size(ysd));
plot(ysdt(i,1),ysdt(i,2),'r');plot(ysdt(end,1),ysdt(end,2),'r*');
plot(ysdt(i,1),ysdt(i,2),'r.');

hold off



display('The same as before with Grand Schmidt') 
% Find conjugate directions

%[x,rho,niter,J,S]=conj_direc(A,b,tol,length(xmin),[0 0]);
theta=pi/2;
rot=[cos(theta) -sin(theta); sin(theta) cos(theta)];
  
u1=b;
u2=rot*u1;

p1=u1/norm(u1);
p2=u2-(u2'*A*p1)/(p1'*A*p1)*p1
p2=p2/norm(p2);

S=[p1 p2];

display('Test for Conjugacy') 
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1


bb=S'*b;
D=S'*A*S;
SI=inv(S);

figure(3)

subplot(221)
contour(x1,x2,z,20); hold on
[ysd,fy,normdeltay,normgrad]=desc_dir(A,b,S,[0 0],maxiter,tol)


i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off

subplot(222)
zz=1/2*(x1.*(D(1,1).*x1+D(1,2).*x2)+x2.*(D(2,1).*x1+D(2,2).*x2))-bb(1).*x1-bb(2).*x2;
contour(x1,x2,zz,20); hold on

ysdt=SI*ysd';
ysdt=ysdt';

i=1:max(size(ysd));
plot(ysdt(i,1),ysdt(i,2),'r');plot(ysdt(end,1),ysdt(end,2),'r*');
plot(ysdt(i,1),ysdt(i,2),'r.');



hold off



return;


% The same as before with conjugate directions preconditioned
% Find conjugate directions

[x,rho,niter,J,S]=conj_direc(A,b,tol,length(xmin),[0 0]);
S=S(:,1:2);
lambda=svd(A);

%p1=1/sqrt(lambda(2))*S(:,1);
%p2=1/sqrt(lambda(1))*S(:,2);

p1=S(:,1);
p2=S(:,2);

p10=p1;
p20=p2;

p1=p1/sqrt(p1'*A*p1);
p2=p2/sqrt(p2'*A*p2);

S=[p1 p2];
SI=inv(S);
%Test for conjugacy and orthgonality
test1=p1'*A*p2
test2=p2'*A*p1
test3=p1'*p2
test4=p2'*p1


bb=S'*b;
D=S'*A*S;

figure(4)

subplot(221)
contour(x1,x2,z,20); hold on
[ysd,fy,normdeltay,normgrad]=desc_dir(A,b,S,[0 0],maxiter,tol)


i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off

subplot(222)
zz=1/2*(x1.*(D(1,1).*x1+D(1,2).*x2)+x2.*(D(2,1).*x1+D(2,2).*x2))-bb(1).*x1-bb(2).*x2;
contour(x1,x2,zz,20); hold on
ysdt=SI*ysd';
ysdt=ysdt';

i=1:max(size(ysd));
plot(ysdt(i,1),ysdt(i,2),'r');plot(ysdt(end,1),ysdt(end,2),'r*');
plot(ysdt(i,1),ysdt(i,2),'r.');

hold off



subplot(223)
[ysd,fy,normdeltay,normgrad]=desc_dir(D,bb,bb(:),[0 0],maxiter,tol)
contour(x1,x2,z,20); hold on
ysdt=SI*ysd';
ysdt=ysdt';

i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off



%i=1:max(size(ysd));
%plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
%plot(ysd(i,1),ysd(i,2),'r.');
%hold off









