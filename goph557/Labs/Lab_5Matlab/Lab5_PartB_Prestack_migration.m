%Lab5_PartB_Prestack_migration
clear
close all

dx=10; % grid size
dz=10;
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
[nz nx]=size(vel);
velsmooth=gaussian_smoother(vel,1:nx,1:nz,10);
veldirectwave=ones(nz,1)*vel(1,:);

figure;
imagesc(xvel,zvel,vel);
title('Velocity Model');
xlabel('Distance (m)');
ylabel('Depth (m)');
colorbar;
colormap('jet');
prepfig;
set(gcf,'position',[100 100 800 400]);

figure;
imagesc(xvel,zvel,veldirectwave);
title('Velocity Model');
xlabel('Distance (m)');
ylabel('Depth (m)');
colorbar;
colormap('jet');
prepfig;
set(gcf,'position',[100 100 800 400]);
%%
x=0:dx:(size(vel,2)-1)*dx;
z=0:dz:(size(vel,1)-1)*dz;

figure;
imagesc(z,x,vel);title('velocity model for comparison with reflectivity');
xlabel('midpoint');ylabel('Depth(m)');prepfig;

snap1=zeros(size(vel));%snap1 is the wavefield before the source goes off
snap2=snap1;
% xmax
% put souce at xmax/3, z=0
% ixs=round(xmax/(3*dx10))+1;
%make a shot record with the 5 point Laplacian
dtstep=.0005;
dt=.004;
tmax=1.5;
fdom=10;%dominant frequency of wavelet
tmaxw=.4;%wavelet length
[w,tw]=ricker(dt,fdom,tmaxw);%make a minimum phase wavelet
laplacian=1;
nshots=2;
dxs=100;

%create shot locations vector;
xshots=zeros(nshots,1);
for is=1:nshots
    ixs=80+is*dxs;
    xshots(is)=ixs*dx;
end

for is=1:nshots
    snap2(5,ixs)=1;%snap2 is zero except at the source location
    xshots=ixs*dx;
    [shotrecord1,seis,t]=afd_shotrec(dx,dtstep,dt,tmax,vel,snap1,snap2,x,zeros(size(x)),w,tw,laplacian);
    [shotrecord_smo,seis,t]=afd_shotrec(dx,dtstep,dt,tmax,veldirectwave,snap1,snap2,x,zeros(size(x)),w,tw,laplacian);
end

figure,imagesc(t,x,shotrecord1);title('First shot to be migrated');
shotrecord2=zeros(size(shotrecord1));
shotrecord2=shotrecord1-shotrecord_smo; % one trace shot record
figure;imagesc(t,x,shotrecord2);title('First shot without direct wave');
%%
%% Prestack Kirchhoff Time Migration
[vrms,tv,vint]=vz2vrms(vel,z,t(2)-t(1),max(t)); % create RMS velocity
[shotmig,tmig,xmig]=kirk_shot(shotrecord2,t,x,xshots,vrms,tv,x);
imagesc(tmig,xmig,shotmig);title('Migration of one shot gather with time Kirchhoff method');
xlabel('midpoint');ylabel('time')
 
%% Prestack PSPI Depth Migration
frange=[1 50];
stab=0.001;
%data=zeros(size(shotrecord1));
%migrate the selected shots
[shotsmigPSPI,shotsmigPSPI_cc,illumination]=pspi_shot(shotrecord2,t,x,vel,x,z,xshots,frange,stab);
figure;
imagesc(z,x,shotsmigPSPI_cc);title('PSPI with CC imaging condition');
xlabel('midpoint');ylabel('Depth(m)');prepfig;
%%