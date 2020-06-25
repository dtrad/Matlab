function [dout,dout2]=filter2d(d,f,F)
% Apply shaping filter
% y: freq domain filtering
% y2: time domain filtering

%close all
%clear all

%load sudata.mat

[n1 n2]=size(d);

D=fft2(d);
%F=fft2(f);

DOUT=D.*F;
DOUT=duplic2d(DOUT);

dout=real(ifft2(DOUT));

% Comment out these two lines if want to filter in freq
dout2=dout;
return



% Time domain 2d deconvolution
dout2=conv2(d,f);
n12=n1/2;
n22=n2/2;

dout2=dout2(n12:n12+n1-1,n22:n22+n2-1);

return;



