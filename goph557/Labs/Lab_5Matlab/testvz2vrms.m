%%
%make an rms velocity model
load matlab
[vrms,tv,vint]=vz2vrms(vel,z,t(2)-t(1),max(t));
%vrms=vint2vrms(vel,t);
%plot

figure;
imagesc(xrec,tv,vint);
title('Interval Velocity in Time ');
xlabel('Distance (m)');
ylabel('Time (s)');
colorbar;
colormap('jet');
prepfig;
set(gcf,'position',[100 100 800 400]);

figure;
imagesc(xrec,tv,vrms);
title('RMS Velocity in Time ');
xlabel('Distance (m)');
ylabel('Time (s)');
colorbar;
colormap('jet');
prepfig;
set(gcf,'position',[100 100 800 400]);


