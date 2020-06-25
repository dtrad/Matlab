function mindescdir(A,b,S,x0,z,x1,x2,nfig,char)
% Given a quadric surface f(x)=1/2 x'Ax + b'x , and a set of search
% directions, this programs :
% 1 - Minimize the quadric, and plot the function and path
% 2 - Rotates the quadric using the search direction axes.
% 3 - Minimize the rotated quadric along the coordinate axes
% 4 - Plot the function and path.
% 5 - Takes the path in the rotated quadric, applies inverse
% rotation to the original system and plot it.
% Input 
%    A Hessian
%    b linear term for the quadric
%    S search vectors arranged in column vector shape ([p1(:) p2(:) ..
%    x0 Staring point
%    z Original quadric (for plotting only) 
%    x1,x2 grid axes (for plotting)
%    nfig : figure number where to plot
%    char legend
% Output
%    
% Daniel Trad - EOSC555 

if (nargin < 9) char='';end 
if (nargin < 8) nfig=1 ;end

locfig=220 % first plot location -1 inside fig nfig

tol=1e-5;  
maxiter=40;

% Use later for inverse transformation
SI=inv(S);

% Rotate quadric
bb=S'*b;
D=S'*A*S;

% Work in original system
figure(nfig)
subplot(locfig+1)
contour(x1,x2,z,20); hold on
[ysd,fy,normdeltay,normgrad]=desc_dir(A,b,S,x0,maxiter,tol)
title('Minimization in the original system')

text2=sprintf('xmin=[%4.2f %4.2f] \n  %d iterations \n',ysd(end,1),ysd(end,2), ...
	      length(ysd)-1);
text(-1.,1.5,text2);

i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off

% Work in rotated system
subplot(locfig+2)
zz=1/2*(x1.*(D(1,1).*x1+D(1,2).*x2)+x2.*(D(2,1).*x1+D(2,2).*x2))-bb(1).*x1-bb(2).*x2;
contour(x1,x2,zz,20); hold on
title('Minimization in the rotated system')
% Before minimization rotate starting point
x0t=SI*x0;
SS=eye(2); % The rotated axes SS=inv(S)*S=idendity
% Minimize in the rotated system along coordinate axes
[ysd,fy,normdeltay,normgrad]=desc_dir(D,bb,SS,x0t,maxiter,tol)
text2=sprintf('xmin=[%4.2f %4.2f] \n  %d iterations \n',ysd(end,1),ysd(end,2), ...
	      length(ysd)-1);
text(-1.,0.5,text2);
% Plot path
i=1:max(size(ysd));
plot(ysd(i,1),ysd(i,2),'r');plot(ysd(end,1),ysd(end,2),'r*');
plot(ysd(i,1),ysd(i,2),'r.');

hold off

% Transform results from rotated quadric to original system
subplot(223)
contour(x1,x2,z,20); hold on
title('Minimization in the x_{rot} system, plotted in the x ')
% Rotate path
ysdt=S*ysd';
ysdt=ysdt';

i=1:max(size(ysd));
plot(ysdt(i,1),ysdt(i,2),'r');plot(ysdt(end,1),ysdt(end,2),'r*');
plot(ysdt(i,1),ysdt(i,2),'r.');
hold off

subplot(221)
text(5,6.5,char)

