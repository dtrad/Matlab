clear;
close all;
clc;

% %% Layered model

% vel_true=ones(100,100)*3000;
% [nz,nx]=size(vel_true);
% vel_true(60:end,:)=4000;
% %% blocky model
% % vel_true(75:100,60:150)=3000;
% 
% density_true=ones(nz,nx)*2000;
% %density_true(51:end,:)=2500;
% 
% dx=8;
% dz=8;
% vel_initial=gaussian_smoother(vel_true,1:nx,1:nz,3);
% vel_initial=ones(nz,nx)*3000;
% density_initial=density_true;
% density_initial=gaussian_smoother(density_true,1:nx,1:nz,4);

%%
%velocity model of marmousi in Madagascar;
load vel_marmousi_Madagascar.mat;
[N M]=size(marmousi_Madagascar);
vel_true=marmousi_Madagascar(1:1:floor(N*1),1:1:floor(M/1));
[nz,nx]=size(vel_true);
clear marmousi_Madagascar;
density_true=1400*ones(nz,nx);
dx=8;
dz=8;
%load bgvel.mat;
%vel_initial=bgvel;
vel_initial=gaussian_smoother(vel_true,1:nx,1:nz,2);
%vel_initial=vel_true;
density_initial=gaussian_smoother(density_true,1:nx,1:nz,1);

%% plot model
zz=0:dz:(nz-1)*dz;
xx=0:dx:(nx-1)*dx;

% zz=zz./1000;
% xx=xx./1000;

vel_true=1./vel_true;
vel_initial=1./vel_initial;

%% parameterization 

%velocity parameterization % v^-1 
model_true(:,1)=vel_true(:);
model_initial(:,1)=vel_initial(:);
model_true(:,2)=1./density_true(:);
model_initial(:,2)=1./density_initial(:);


%% plot model

imagesc(xx,zz,1./reshape(model_true(:,1),nz,nx));
shading interp
title('True model','Fontname','Consulas','fontsize',18,'fontweight','bold');
xlabel('Distance (m)','Fontname','Consulas','fontsize',18);
ylabel('Depth (m)','Fontname','Consulas','fontsize',18);
colorbar;
axis ij;

imagesc(xx,zz,1./reshape(model_initial(:,1),nz,nx));
shading interp
title('Initial model','Fontname','Consulas','fontsize',18,'fontweight','bold');
xlabel('Distance (m)','Fontname','Consulas','fontsize',18);
ylabel('Depth (m)','Fontname','Consulas','fontsize',18);
caxis([min(1./model_true(:,1)),max(1./model_true(:,1))])
colorbar;
axis ij;


%% Acquisition setting
sz=2; % souce depth
rz=2; % receiver depth

s_inter=100; % source interval
s_initial=5; % initial source position
s_end=nx-2; % end source position
xs=s_initial:s_inter:s_end; % x-coordinates of the sources


r_inter=3; % receiver interval
r_initial=1; % initial receiver position
r_end=nx; % end receiver position
xr=r_initial:r_inter:r_end; % x-coordinates of the receivers

%% Ricker Wavelet
dt=0.001 ;
maxt=3;
t=0:dt:maxt;
t=-maxt/2:dt:maxt/2;
fdom=15; % dominant frequency
[wavelet,tw_ricker]=ricker(dt,fdom,maxt);
wavelet_f=fft(wavelet);
nt=size(t,2);
df=1./((nt-1)*dt)
nyquist=1./(2*dt);
f  = 0:df:nyquist; 
nf = length(f);


%figure,plot(f,abs(wavelet_f(1:length(f))))
%title('Ricker wavelet');
%xlabel('f(Hz)','fontsize',20);
%ylabel('A','fontsize',20);
%prepfig;

%% minimum phase wavelet
% dt=0.001;
% maxt=1;
% t=0:dt:maxt;
% fdom=15; 
% [wavelet_mini,tw_mini]=wavemin(dt,fdom,maxt);
% wavelet_f=fft(wavelet_mini);
% nt=size(t,2);
% f  = 0:1/((nt-1)*dt):.5/dt; 
% nf = length(f);
% df=f(2)-f(1);
% 
% 
% plot(f,abs(wavelet_f(1:nf)))
% title('minimum phase wavelet');
% axis([0,100,0,1.2]);
% xlabel('f(Hz)','fontsize',20);
% ylabel('A','fontsize',20);
% prepfig;

% plot(t,wavelet_mini);
% title('minimum phase wavelet');
% xlabel('t(s)','fontsize',20);
% ylabel('A','fontsize',20);
% prepfig;
RTM2T=zeros(size(vel_initial));
RTMT=zeros(size(vel_initial));
LSRTMT=zeros(size(vel_initial));
%% Frequency for inversion
maxfreq=20;
imaxf=20/df;
DATAT=zeros(length(xr),length(xs),imaxf);

% for ifreq=2:imaxf
% % frequency=[f_index:0.2:f_index+5];
%     frequency=f(ifreq:ifreq);
%     printf('frequency = %d\n',frequency);    
%     Data = FDFD2D_vd_parameterization_encoding(model_true,frequency,nz,nx,dz,dx,xs,xr,sz,rz,wavelet_f,f);
%     DATAT(:,:,ifreq)=Data;
% end
% save data.mat DATAT xr xs



for ifreq=2:imaxf
% frequency=[f_index:0.2:f_index+5];
frequency=f(ifreq:ifreq);
fprintf('# calculating ifreq %d\n',frequency);

%% Forward modelling

%F = @(model_true)FDFD2D_vd_parameterization_encoding(model_true,frequency,nz,nx,dz,dx,xs,xr,sz,rz,wavelet_f,f);
Data = FDFD2D_vd_parameterization_encoding(model_true,frequency,nz,nx,dz,dx,xs,xr,sz,rz,wavelet_f,f);
DATAT(:,:,ifreq)=Data;

%% Inversion Hessian vector product method
cgmax=5;
Inversion = @(x)Inversion_vd_parameterization_encoding_pwy(x,Data,frequency,nz,nx,dz,dx,xs,xr,sz,rz,wavelet_f,f);
[RTM,LSRTM,RTM2]=Optimization_LSRTM(Inversion,model_initial,model_true,nz,nx,cgmax);
RTMT=RTMT+reshape(RTM,nz,nx);
LSRTMT=LSRTMT+reshape(LSRTM,nz,nx);
RTM2T=RTM2T+reshape(RTM2,nz,nx);
end
save data.mat DATAT xr xs
RTM=RTMT;
LSRTM=LSRTMT;
RTM2=RTM2T;
clear RTMT LSRTMT RTM2T
%% plot result

RTM=reshape(RTM,nz,nx);
LSRTM=reshape(LSRTM,nz,nx);
%RTM=RTM./vel_initial;
%LSRTM=LSRTM./vel_initial;
plot_result1(xx,zz,RTM,LSRTM,nz,nx);

figure,imagesc(RTM);
figure,imagesc(LSRTM);

writeBin('rtm2',RTM2);
writeBin('rtm',RTM);
writeBin('lsrtm',LSRTM);
writeBin('velocity',1./vel_initial);
% imagesc(xx,zz,-reshape(RTM,nz,nx));
% title('RTM image','Fontname','Times New Roman','fontsize',18);
% xlabel('Distance (km)','Fontname','Times New Roman','fontsize',18);
% ylabel('Depth (km)','Fontname','Times New Roman','fontsize',18);
% colormap(gray);
% 
% imagesc(xx,zz,reshape(LSRTM,nz,nx));
% title('LSRTM image','Fontname','Times New Roman','fontsize',18);
% xlabel('Distance (km)','Fontname','Times New Roman','fontsize',18);
% ylabel('Depth (km)','Fontname','Times New Roman','fontsize',18);


%% (laplacian) filter
% % 
%   h1=[0 1 0;1 -4 1;0 1 0];
%   h2=[0 1 0;1 -4 1;0 1 0];
%  RTM_f=imfilter(reshape(RTM,nz,nx),h1,'replicate');
%  LSRTM_f=imfilter(reshape(LSRTM,nz,nx),h2,'replicate');
% %   h3=[0 2 0;0 0 0;0 -2 0];
% %  LSRTM_f=imfilter(reshape(LSRTM,nz,nx),h3,'replicate');
%  
%  clip1=[-1e13 1e13];
%  clip2=[-5e-5 5e-5];
% % plot_result1(xx,zz,RTM_f,LSRTM_f,nz,nx);
%  plot_result(xx,zz,RTM_f,LSRTM_f,nz,nx,clip1,clip2);





