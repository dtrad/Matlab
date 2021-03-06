function [d2,d,d0,wd,twd]=monsterpiece(perc,method,noise);
echo on
% Example of denoising with 2D wavelet transform and ridgelets
% usage [d2,d,d0,wd,twd]=slantline(perc,method);
% Perc is the number of coefficients to keep if > 1
% or the percentage is < 1.
% It needs the data in file slanttest2 that is a simple 
% line. An ps figure is created for both examples.
% ~/ps/wavedenoise.ps
% ~/ps/ridgedenoise.ps

if (nargin < 3) noise=0.;end
if (nargin < 2) method='ridgelet';end
if (nargin < 1) perc=0.99;end


clip=[-0.0001 0.0001];
clip2=[-.1 .1];

close all 
method=method(1:7);
fig=0

N=256
[x]=readsudata('/home/dtrad/work/datagroll.bin',N,N);

%x=real(x');

d=real(x)/max(max(abs(x)));;
[nt nh]=size(d);

%d=cumsum(d);
d0=d;
d=d+noise*randn(size(d));

figure(1); 
subplot(321); imagesc(d0,clip2);colorbar;
title('data');
title('(a)');

subplot(322),imagesc(d,clip2);colorbar;
title('noisy data');
title('(b)');

if (method=='wavelet')
    percw=perc;
    [d2,wd,twd]=applywt2(d,nt,nh,perc);
    save wavelet_result3.mat d2 wd twd d d0 percw;

    subplot(323); imagesc(wd);colorbar;title('wavelet coefficient before thresholing');
    subplot(324); imagesc(twd);colorbar;title('wavelet coefficients after theresholing');
    text1=sprintf('wavelet denoised data (%d coeff)',perc)
    subplot(325); imagesc(d2);colorbar;title(text1);
    if (fig) print -dps2 ../ps/wavedenoise.ps,end
elseif (method=='ridgele')
    percr=perc;


    mask=zeros(2*nt,2*nh);
    mask(:,nh:end)=1;
    [d2,wd,twd]=applyridget(d,perc,mask);
    d2r=d2;wdr=wd;twdr=wdr;
    
    save ridge_result3.mat d2r wdr twdr d d0 percr ;

    subplot(323); imagesc(wd,clip);colorbar;
    %title('ridgelet coeff. before thresholding');
    title('(c)');

    subplot(324); imagesc(twd,clip);colorbar;
%   title('ridgelet coeff. after thresholding');
    title('(d)');
    
    %ext1=sprintf('ridgelet denoised data (%d coeff)',perc)
    subplot(325); imagesc(d2,clip2);colorbar;
    %itle(text1);
    title('(e)');
    
    subplot(326); imagesc(d-d2,clip2);colorbar;
    %itle('cleaned data.');
    title('(f)');
    if (fig) 
      print -dps2 ../ps/ridgedenoise.ps,
    end
    
    if (1) 
      figure,
      subplot(221);mesh(d);title('(a)');VV=axis;
      subplot(222);mesh(d2);title('(b)');axis(VV);
      subplot(223);mesh(d-d2);title('(c)');axis(VV);
      
    end
    
    
    % If combined figure is required for both method set example
elseif (method=='figures')
    load wavelet_result;
    load ridge_result;
    figure(2);
    v=[-8e-1 8e-1];
    subplot(221); imagesc(d0);colorbar;title('data');
    subplot(222),imagesc(d);colorbar;title('noisy data');
    text1=sprintf('wavelet denoised data (%d coeff)',percw)
    subplot(223); imagesc(d2);caxis(v);colorbar;title(text1);
    text1=sprintf('ridgelet denoised data (%d coeff)',percr)
    subplot(224); imagesc(d2r);caxis(v);colorbar;title(text1);
else
    display('wrong option');
    method
    display('it should be one of wavelet or ridgelet');

end

% write data for su 
if (1) 
  writesudata('/home/dtrad/work/datagrollrid.b.bin',d);
  writesudata('/home/dtrad/work/datagrollrid.noise.bin',d2);
  writesudata('/home/dtrad/work/datagrollrid.rec.bin',d-d2);
  writesudata('/home/dtrad/work/datagrollrid.ridge.bin',wd);
  writesudata('/home/dtrad/work/datagrollrid.ridgethresh.bin',twd);
end

return





