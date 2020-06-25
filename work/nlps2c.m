clear;
close all

x=-2:.2:2;
y=-1:.2:3;
[xx,yy]=meshgrid(x,y);
zz=100*(yy-xx.^2).^2+(1-xx).^2;
figure, mesh(x,y,zz)

% Set up the appropriate colormap
% In this case, the colormap has been chosen to give the surf plot
% a nice healthy banana color.
hsv2=hsv;
hsv3=[hsv2(11:64,:); hsv2(1:10,:)];
% draw the surf plot
surfHndl=surface(x,y,zz,'EdgeColor',[.8 .8 .8]);
axis off;
view(10,55);
colormap(hsv3);
hold on;
[c,contHndl]=contour3(x,y,zz+50,[100 500],'k');
set(contHndl,'Color',[.8 .8 .8]);
drawnow
plot3(-1.9,2,267.62,'ko', ...
      'MarkerSize',15, ...
      'LineWidth',2, ...
      'EraseMode','none');
text(-1.9,2.2,267.62,'   Begin', ...
     'Color',[0 0 0], ...
     'EraseMode','none');
plot3(1,1,0,'ko', ...
      'MarkerSize',15, ...
      'LineWidth',2, ...
      'EraseMode','none');
text(0.8,1.4,0,'   End', ...
     'Color',[0 0 0], ...
     'EraseMode','none');
x=[-1.9 2];

str=[' Steepest Descent'];


%x=[4*rand-2 4*rand-1];
OPTIONS=0;
OPTIONS(1)=-1;
OPTIONS(6)=2;
OPTIONS(10)=200;
OPTIONS(14)=300;
GRAD='[100*(4*x(1)^3-4*x(1)*x(2))+2*x(1)-2; 100*(2*x(2)-2*x(1)^2); xpbanplt(x)]';
f='100*(x(2)-x(1)^2)^2+(1-x(1))^2';
str2=' [x, options] = fminu(f,x,OPTIONS,GRAD);';
[x, options] = fminu(f,x,OPTIONS,GRAD);