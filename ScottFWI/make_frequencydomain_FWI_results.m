clear;
clc;
close all;

%% create velocity model
tic;

Ball_model;
%Layer_model;

% Thickness of PML layer, and theoretical reflection coefficient from model
% boundary. If R_theor goes down, PML_thick should go up to prevent
% reflections from numerical errors
PML_thick = 9;
R_theor = 1e-3;
% PML_thick = 21;
% R_theor = 1e-4;

NPML = (nz+2*PML_thick)*(nx+2*PML_thick);

% Velocities are input as vectors
vel0 = vel_initial(:);
vel_true = vel_true(:);
vmin = min(vel_true);
vmax = max(vel_true);
%% Acquisition

% Source and receiver coordinates should be in grid number
s_inter=1; % source interval
s_initial=2; % initial source position
s_end=nx-2; % end source position
sx=s_initial:s_inter:s_end; % x-coordinates of the sources
sz = 3*ones(size(sx));

% Example of sources at top and bottom
%sz = [3*ones(size(sx)), 0.5*nz*ones(size(sx)), (nz-2)*ones(size(sx))];
%sx = [sx,sx,sx];
sz = [3*ones(size(sx)),  (nz-2)*ones(size(sx))];
sx = [sx,sx];
ns = length(sx);
ns

r_inter=1; % receiver interval
r_initial=2; % initial receiver position
r_end=nx-1; % end receiver position
rx=r_initial:r_inter:r_end; % x-coordinates of the receivers
rz = 2*ones(size(rx));

[S,R] = Define_Acquisition(sz, sx, rz, rx, nz, nx, PML_thick);

figure(9998);
imagesc(reshape(vel_true,nz,nx));
hold on;
plot(sx,sz,'ko');
plot(rx,rz,'k+');
caxis([vmin,vmax]);
title('True model')
drawnow

%% Determine frequencies used
% Number of frequency bands inverted
numbands = 13;  % 7;   %13;
finc = 2;
% Number of frequencies per band
step=5;
freq=zeros(1,numbands*step);
%startfreq = ones(1,numbands);
%endfreq = linspace(1,25,numbands);
f1 = 1;
f9 = floor(f1 + (numbands-1)*finc);  % 4);    %*2);
startfreq = f1*ones(1,numbands);
endfreq = linspace(f1,f9,numbands);    %5 + (13-1)*2
%endfreq = linspace(5,25,numbands);

for n=1:numbands
    freq(1 + (n-1)*step : n*step) = linspace(startfreq(n), endfreq(n), step);
    fprintf('n %d freq %f to %f\n',n,freq(1+(n-1)*step),freq(n*step));
end

% Amplitude of wavelet at each frequency in freq
fwave=ones(1,length(freq));

% fwave = exp(-((freq - 20)/10).^2);

%% Finite difference function
FDFDfunc = @(frequency,fwave,vel)FDFD(vel,vel0,frequency,S,fwave,PML_thick,R_theor,nz,nx,dz,dx);

%% Generate data
D = Get_data( FDFDfunc, freq, fwave, vel_true, R );

%% Run inversion
% Number of iterations per band
numits=1;
% Optimization used: 1 - Steepest Descent, 2 - Gauss Newton
optype = 2;
% 
rangevel = 1.1*max(abs(vel_true-vel0));
% Main FWI code
vel = FDFWI( D, freq, step, fwave, FDFDfunc, nz, nx, vel0, R, optype, numits, PML_thick, rangevel,vmin,vmax,vel_true,sx,sz,rx,rz);

%{
figure(9998);
imagesc(reshape(vel_true,nz,nx));
caxis([vmin,vmax]);
title('True model')
drawnow
%}

%% Plot profiles
% close all;

profileindex = floor(nx*0.5);  %25;
velmat = reshape(vel,nz,nx);
velprofile = velmat(:,profileindex);
veltruemat = reshape(vel_true,nz,nx);
veltrueprofile = veltruemat(:,profileindex);
          
t_cur = toc;
t_h = floor(t_cur/3600);
t_m = floor((t_cur-3600*t_h)/60);
t_s = round(t_cur - 3600*t_h - 60*t_m);
if (optype==2)
 txt = sprintf('Gauss-Newton: Freq: %f to %f by %f\n Elapsed time: %i minutes and %i seconds',f1,endfreq(numbands),finc,t_m,t_s);
else
 txt = sprintf('Steepest Descend: Freq: %f to %f by %f\n Elapsed time: %i minutes and %i seconds',f1,endfreq(numbands),finc,t_m,t_s);  
end    
txt_time=sprintf('Elapsed time: %i minutes and %i seconds',t_m,t_s);

%close(9999);
figure(9999);
subplot(2,2,[1 2]), 
plot(1:nz,velprofile,'k-',1:nz,veltrueprofile,'b-');
title(txt);
ylim([1800 3400]);
grid; grid minor;
%txt=sprintf('Freq: %f to %f by %f',f1,endfreq(numbands),finc);

subplot(2,2,3),
imagesc(velmat); %colormap('gray'); 
caxis([vmin,vmax]);
colorbar;

subplot(2,2,4),
imagesc(veltruemat); %colormap('gray'); 
hold on;
plot(sx,sz,'ko');
plot(rx,rz,'k+');
caxis([vmin,vmax]);
colorbar;


    