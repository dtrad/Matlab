%
% Takes the adjoint obtained by HRT and the desired sparse RT, 
% calculates a shaping filter to transform adjoint into sparse RT
% and applies this filter to other events.
%
% Daniel Trad - CDSST - UBC August, 2002

close all
clear all
% Read the data
% load sudata.mat

[d]=readsudata('/home/dtrad/work/adjoint.bin',512,100);
[dout]=readsudata('/home/dtrad/work/modells.bin',512,100);

% d adoint 
% dout desired sparse RT

% calculate filter using a window
f1=185
n1=64
f2=1
n2=100


dw=window(d,f1,n1,f2,n2);
doutw=window(dout,f1,n1,f2,n2);

[f,F,y,y2]=shapingfilter2d(dw,doutw,.0001);

figure(1)
subplot(221);simage(doutw);colorbar;title('desired output')
subplot(222);simage(dw);colorbar;title('input')
subplot(223);simage(f);colorbar;title('filter')
subplot(224);simage(y);colorbar;title('filtered input FD')

figure(2)
subplot(221);simage(doutw);colorbar;title('desired output')
subplot(222);simage(dw);colorbar;title('input')
subplot(223);simage(f);colorbar;title('filter')
subplot(224);simage(y2);colorbar;title('filtered input TD')

% Apply the filter in windows in the data set
[dz]=filter2dwin(d,f,F);

figure(3)
subplot(221);simage(d);title('Adjoint');
subplot(222);simage(dout);title('LS HRT');
subplot(223);simage(dz);title('Filtered Adjoint');
subplot(224);simage(dout-dz);title('difference');


writesudata('/home/dtrad/work/modelfilt.bin',dz)

