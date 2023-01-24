function plot_result(xx,zz,RTM,LSRTM,nz,nx,clip1,clip2)


figure;
subplot(1,2,1)
pcolor(xx,zz,-RTM);
title('RTM image','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (km)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (km)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;
caxis(clip1);
colormap(gray);
subplot(1,2,2)
pcolor(xx,zz,LSRTM);
title('LSRTM image','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (km)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (km)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;
caxis(clip2);
colormap(gray);
colormap(rb);
prepfig;