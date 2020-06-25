function [yd,y,wy,twy]=ridgegroll(s,perc)
% Ground roll prediction by thresholding in the ridgelet domain
% 
% ridgegroll(filename,perc)
%
% reads file "filename.bin" and writes the file "filename.noise.bin"
% in the path given by 'mypath'
%
% It needs BeamLab and writesudata.m
%
% For groundroll removal is very important that the left hand 
% side of the ridgelet plane is zeroed because signal maps mainly 
% to this half plane (basically horizontal maps left, basically
% vertical maps right.
%
% D. Trad and Felix Herrmann - Nov, 2002

mypath='/home/dtrad/work/'

clip=[-0.2,0.2];

close all
t0=129;N=128;
if (nargin<1) perc=0.005;end
% load data

[d]=readsudata([mypath,s,'.bin'],128,86);
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
imagesc(y,clip);title('data');

% Apply ridgelet transform
wy=FastOrthoRidgeletTransform(y);
save wy.mat wy

%load wy.mat
figure(2);imagesc(wy);title('ridgelet coefficients');colorbar;
figure(gcf); 

% Thresholding 
q=threshold(wy,perc);
twy=wy;
twy(find(abs(wy)<q))=0;
twy(:,1:128)=0;
% Use a mask from a model of the noise
%load ttwy.mat;
%twy=wy.*ttwy;

%ttwy=zeros(size(twy));
%ttwy(1:64,1:64)=twy(1:64,1:64);


figure(3);imagesc(twy);title('thresholded ridgelet coefficients');colorbar;
figure(gcf)
% Inverse ridgelet transform
yn=Inv_FastOrthoRidgeletTrans(twy);
yn=real(yn);

figure(4);imagesc(yn,clip);title('Noise');
figure(gcf)

yd=y-yn;
figure(5);imagesc(yd,clip);title('Denoised data');
figure(gcf)


writesudata([mypath,s,'.noise.bin'],yn);

return




