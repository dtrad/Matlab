
subplot(2,2,1);
imagesc(xx,zz,-RTM);
title('RTM image FD','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (m)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (m)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;

subplot(2,2,2)
imagesc(xx,zz,LSRTM);
title('LSRTM image FD','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (m)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (m)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;

subplot(2,2,3)
imagesc(xx,zz,rtmf);
title('RTM image TD','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (m)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (m)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;

subplot(2,2,4)
imagesc(xx,zz,lsrtmf);
title('LSRTM image TD','Fontname','Times New Roman','fontsize',18);
xlabel('Distance (m)','Fontname','Times New Roman','fontsize',18);
ylabel('Depth (m)','Fontname','Times New Roman','fontsize',18);
shading interp
axis ij;
colormap('gray')


