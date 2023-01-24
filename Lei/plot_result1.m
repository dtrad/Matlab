function plot_result1(xx,zz,RTM,LSRTM,nz,nx)


figure;
subplot(1,2,1)
pcolor(xx,zz,-RTM);
title('RTM image ','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (m)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (m)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;
colormap(gray)
subplot(1,2,2)
pcolor(xx,zz,LSRTM);
title('LSRTM image','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (m)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (m)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;
colormap(gray);
%prepfig;