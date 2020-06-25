%function [yd,y]=ridgepacketsgroll(D,perc)
clear all
close all
D=[7 7];
perc=100;
% Denoise by thresholding in the ridgelet domain
% testridgeden(perc,noise)
close all
%t0=129;N=256;
%if (nargin<1) perc=0.5;end
% load data
%load ozdata
%y=[d d d d];
%y=y(t0:2:t0+2*N-1,1:N);

load lenna;

y=lenna(1:128,1:128);
[N M]=size(y);
figure(1)
imagesc(y);title('data');

% Apply ridgelet packet transform

RPFT= fft2_rp(y);
ARPkt = CalcRPPktTable(RPFT,D,'Sine');
RPTree = CalcRPStatTree(ARPkt,D,'Entropy');




[btree] = BestRPBasis(RPTree,D);


coef  = FPT2_RPkt(btree,y,D);
figure(2);imagesc(coef);title('ridgelet coefficients');colorbar;
figure(gcf); 

img0 = IPT2_RPkt(btree,coef,D);
figure;imagesc(img0);title('img0 (all coefficients');

if (0)
    % Thresholding 
    % Define the threshold such that only 100 coefficients remain
    [N1 M1]=size(coef);
    total=N1*M1;
    kept=round((total)*(perc/100));
    if (kept>total) 
        kept=total;
    end
    message=sprintf('number of coefficients to keep=%d from %d or %f percent',...
        kept,total,kept/total*100)



    ss=sort((abs(coef(:))));q=ss(length(ss)-kept+1);
    tcoef=coef;tcoef(find(abs(coef)<q))=0;
    figure(3);imagesc(tcoef);title('thresholded ridgelet coefficients');colorbar;
    figure(gcf)
    img = IPT2_RPkt(btree,tcoef,D);
else
    img=img0;
end




figure;imagesc(y);title('y');
figure;imagesc(img);title('img')
figure;imagesc(y-img);title('y-img');

figure,imagesc(log(fftshift(abs(fft2(y)))));colorbar;title('y');
figure,imagesc(log(fftshift(abs(fft2(img)))));colorbar;title('img');
figure,imagesc(log(fftshift(abs(fft2(y-img)))));colorbar;title('y-img');



return;

