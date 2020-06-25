function [yd,y] = curvegroll(filename2,perc,pfilt,dfilt,dlevs,rlevs,par)
% Ground roll prediction by thresholding in the ridgelet domain
% reads file "filename.bin" and writes the file "filename.noise.bin"
% in the path given by 'mypath'
% 
% curvegroll(filename,perc,pfilt,dfilt,rlevs,par)
%
% pfilt and dfilt define the wavelet
% dlevs define the levels to use in the contourlet decomposition
% Examples
% pfilt     =   'db8';
% dfilt     =   'pkva';
% dlevs     =   [1,1,2,2];
% rlevs is a cell array that defines a mask. Only those values of
% rlevs=1  will serve in  the prediction.  
% Example 
% rlevs           =   cell(length(dlevs),1);
% rlevs{5}{3}     =   0;
%
% It needs BeamLab, plotting routines in GT2.9a and writesudata.m
%
% For groundroll removal is very important that the some quadrants 
% of the countourlet plane where the signal maps are zeroed because
% mainly signal maps to these quadrants
%
% Felix Herrmann and D. Trad - Nov, 2002
%


clip=[-0.2,0.2];
mypath='/home/dtrad/work/'


close all
t0=129;N=128;
if (nargin<1) perc=0.005;end
% load data

[mypath,filename2,'.bin']
[d]=readsudata([mypath,filename2,'.bin'],128,86);

%d = su2mat(s,[128,86]);

dd=[d zeros(size(d))];
d=dd(:,1:128);
%d=d.';
%load ozdata1win

%y=y(t0:2:t0+2*N-1,1:N);

% Apply a low pass filter

if (0) 
  [B,A]=butter(2,0.01);
  dd=zeros(size(d));
  for ih=1:N
    dd(:,ih)=filtfilt(B,A,d(:,ih));
  end
  y=dd;
  
else
  y=d;
end


[N M]=size(y)
figure(1)
imagesc(y,clip);colormap(seiscolor);title('data');

% Apply curvelet transform


%wy=FastOrthoRidgeletTransform(y);
wy =  pdfbdec(y,pfilt,dfilt,dlevs);
save wy.mat wy

%load wy.mat

figure(2);
[h] = ImageCurveletmod(wy,dlevs,1);
axes(h(1));
title('curvelet coefficients');
figure(gcf); 

% $$$ % Thresholding 
% $$$ [vwy, s]              = pdfb2vec(wy);
% $$$ q        = threshold(vwy,perc);
% $$$ 
% $$$ tvwy                  = vwy;
% $$$ tvwy(find(abs(vwy)<q))= 0;

% set the subbands to zero

if nargin > 6,
  for istr=1:length(dlevs)
    dstring(istr)=num2str(dlevs(istr));
  end
    
  if ~exist([mypath,'nvar_' dstring '.mat'],'file')
    nvar = pdfb_nest(N, pfilt, dfilt, dlevs);
    save([mypath,'nvar_' dstring '.mat'],'nvar');
  else
    nvar =  load([mypath,'nvar_' dstring '.mat']);
    
  end

  twy                   =  ThresBand(wy,dlevs,rlevs,par,size(d,1),nvar.nvar);
else
  twy                   =  ThresBand(wy,dlevs,rlevs);
end
%keyboard
figure(3);
[h]=ImageCurvelet(twy,dlevs,1);
axes(h(1));
title('thresholded curvelet coefficients');
figure(gcf)

% Inverse curvelet transform

yn   =  pdfbrec(twy  ,pfilt, dfilt);

figure(4);imagesc(yn,clip);colormap(seiscolor);title('Noise');
figure(gcf)

yd=y-yn;
figure(5);imagesc(yd,clip);colormap(seiscolor);title('Denoised data');
figure(gcf)


% save the noise for later subtraction
writesudata([mypath,filename2,'.noise.curve.bin'],yn);

return




