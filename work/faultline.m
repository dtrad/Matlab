function [d2,d]=faultline(perc,method);
% Example of denoising with 2D wavelet transform and ridgelets
close all 
method=method(1:7);
fig=1

load faultline
[nt nh]=size(d);
d0=d;

d=d+randn(size(d));

figure(1); 
subplot(321); imagesc(d0);colorbar;title('data');
subplot(322),imagesc(d);colorbar;title('noisy data');

if (method=='wavelet')
    [d2,wd,twd]=applywt2(d,nt,nh,perc);
    save wavelet_result.mat d2 wd twd d d0;

    subplot(323); imagesc(wd);colorbar;title('wavelet coefficient before thresholing');
    subplot(324); imagesc(twd);colorbar;title('wavelet coefficients after theresholing');
    text1=sprintf('wavelet denoised data (%d coeff)',perc)
    subplot(325); imagesc(d2);colorbar;title(text1);
    if (fig) print -dps2 wavedenoise.ps,end
elseif (method=='ridgele')
    [d2,wd,twd]=applyridget(d,nt,nh,perc);
    d2r=d2;wdr=wd;twdr=wdr;
    
    save ridge_result.mat d2r wdr twdr d d0;

    subplot(323); imagesc(wd);colorbar;title('ridgelet coeff. before thresholding');
    subplot(324); imagesc(twd);colorbar;title('ridgelet coeff. after thresholding');
    text1=sprintf('ridgelet denoised data (%d coeff)',perc)
    subplot(325); imagesc(d2);colorbar;title(text1);

    if (fig) print -dps2 ridgedenoise.ps,end
else
    display('wrong option');
    method
    display('it should be one of wavelet or ridgelet');

end

return
