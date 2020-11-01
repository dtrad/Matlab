%kirchhoff time migration and PSPI depth migration
clear
close all

dt=0.002;
tmax=2;
t=0:dt:tmax;

dtstep=.001; % size of time step for modelling (in seconds)
dx=5; % grid size
dz=5;
velmodel=1


if (velmodel==1)
    [vel,xvel,zvel]=thrustmodel(dx);
elseif (velmodel==2)
    [vel,xvel,zvel]=channelmodel(dx);
elseif (velmodel==3)
    [vel,xvel,zvel]=wedgemodel(dx);
elseif (velmodel==4)
    [vel,xvel,zvel]=flatmodel(dx);
elseif (velmodel==5)
    [vel,xvel,zvel]=synclinemodel(dx);      
end

figure;
imagesc(xvel,zvel,vel);
title('Velocity Model');
xlabel('Distance (m)');
ylabel('Depth (m)');
colorbar;
colormap('jet');
prepfig;
set(gcf,'position',[100 100 800 400]);

%%   
xrec=0:dx:(size(vel,2)-1)*dx;
z=0:dx:(size(vel,1)-1)*dx;
zrec=zeros(size(xrec)); % z-positions of receivers (in consisent units)
filt=[5 10 40 50];
laplacian=2;
boundary=2;
phase=0;

fdom=10;
[w,tw]=wavemin(dt,fdom,tmax);

[seismogram,seis,t]=afd_explode(dx,dtstep,dt,tmax,vel,xrec,zrec,w,tw,laplacian,boundary);

figure;
imagesc(xrec,t,seismogram);
title('ZOS in Time ');
xlabel('Distance (m)');
ylabel('Time (s)');
colormap('gray');
prepfig;
set(gcf,'position',[100 100 800 400]);




%%
%make an rms velocity model
%load matlab (for testing without runnning the FD
[vrms,tv,vint]=vz2vrms(vel,z,t(2)-t(1),max(t));

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


%% Kirchhoff Time Migration 
params=nan*ones(1,12);
[seismig,tmig,xmig]=kirk_mig(seismogram,vrms,t,xrec,params);

figure;
imagesc(xmig,tmig,seismig);
title('Kirchhoff Time Migration Result in Time ');
xlabel('Distance (m)');
ylabel('Time (s)');
colormap('gray');
prepfig;
set(gcf,'position',[100 100 800 400]);

%%

[seismigz,z1]=stretch(seismig,tmig,vel,z,dx,1);
n=find(z1==2500);
figure;
imagesc(xmig,z1(1:n),seismigz(1:n,:));
title('Kirchhoff Time Migration Result in Depth ');
xlabel('Distance (m)');
ylabel('Depth (m)');
colormap('gray');
prepfig;
set(gcf,'position',[100 100 800 400]);


%%
[seismig_z]=kirk_migz(seismogram,vint,dt,dx,dz,params);
frange=[0 100];
zcheck=100:200:2500;

figure;
imagesc(xrec,z,seismig_z);
title('Kirchhoff Depth Migration Result');
xlabel('Distance (m)');
ylabel('Depth (m)');
colormap('gray');
prepfig;
set(gcf,'position',[100 100 800 400]);

%%
%pspi depth migration
pspimig=pspi_stack(seismogram,t,xrec,vel,xvel,zvel,[0 80]);
imagesc(xvel,zvel,pspimig);
title('PSPI depth migration (using pspi_stack)')

