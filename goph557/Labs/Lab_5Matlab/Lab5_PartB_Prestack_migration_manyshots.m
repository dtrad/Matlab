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
nshots=3;
dxs=50;

%create shot locations vector;
ixshots=zeros(nshots,1);
for is=1:nshots
    ixs=80+is*dxs;
    ixshots(is)=ixs;
end
xshots=ixshots*dx;
nt=round(tmax/dt+1);
data=zeros(nt,nx,nshots);
for is=1:nshots
    snap1=zeros(size(vel));snap2=snap1;
    snap2(5,ixshots(is))=1;%snap2 is zero except at the source location        
    [shotrecord,seis,t]=afd_shotrec(dx,dtstep,dt,tmax,vel,snap1,snap2,x,zeros(size(x)),w,tw,laplacian);
    snap1=zeros(size(vel));snap2=snap1;
    snap2(5,ixshots(is))=1;    
    [directwave,seis,t]=afd_shotrec(dx,dtstep,dt,tmax,veldirectwave,snap1,snap2,x,zeros(size(x)),w,tw,laplacian);
    data(:,:,is)=(shotrecord-directwave);
    figure,imagesc(t,x,shotrecord);title('Full shot to be migrated');
    figure,imagesc(t,x,directwave);title('direct wave to be migrated');
    figure,imagesc(t,x,data(:,:,is));title('One shot to be migrated');    
end

%%
%% Prestack Kirchhoff Time Migration
[vrms,tv,vint]=vz2vrms(vel,z,t(2)-t(1),max(t)); % create RMS velocity
migkir=zeros(nt,2*nx-1);
for is=1:nshots
    [shotmig,tmig,xmig]=kirk_shot(data(:,:,is),t,x,xshots(is),vrms,tv,x);
    figure;imagesc(tmig,xmig,shotmig);title('Migration of one shot gather with time Kirchhoff method');
    xlabel('midpoint');ylabel('time')
    migkir=migkir+shotmig;
end
figure,imagesc(tmig,xmig,migkir);
title('Migration of all shots with time Kirchhoff method');
xlabel('midpoint');ylabel('time');

%% Prestack PSPI Depth Migration
frange=[1 50];
stab=0.001;
migPSPI=zeros(size(vel));
%data=zeros(size(shotrecord1));
%migrate the selected shots
for is=1:nshots
    [shotmigPSPI,shotmigPSPI_cc,illumination]=pspi_shot(data(:,:,is),t,x,vel,x,z,xshots(is),frange,stab);
    figure;
    imagesc(z,x,shotmigPSPI_cc);title('One shot PSPI with CC imaging condition');
    xlabel('midpoint');ylabel('Depth(m)');prepfig;
    migPSPI=migPSPI+shotmigPSPI_cc;
end
figure;imagesc(z,x,migPSPI);
title('Migration of all shots with time Kirchhoff method');
xlabel('midpoint');ylabel('time');
%%