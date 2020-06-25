%function [yd,y]=ridgepacketsgroll(D,perc)
echo on
clear all
close all
warning off

D=[6 6];

% Denoise by thresholding in the ridgelet domain
% testridgeden(perc,noise)
close all
t0=129;N=64;
%if (nargin<1) perc=0.5;end
% load data
load ozdata
y=[d d d d];
y=y(t0:1:t0+N-1,1:N);


[N M]=size(y)
figure(1)
imagesc(y);title('data');

% Apply ridgelet packet transform

RPFT= fft2_rp(y);

ARPkt = CalcRPPktTable(RPFT,D,'Sine');

%load ARPkt.mat

%return
RPTree = CalcRPStatTree(ARPkt,D,'Entropy');


figure(2);imagesc(RPTree);title('ridgelet coefficients');colorbar;
figure(gcf); 


[btree] = BestRPBasis(RPTree,D);


coef  = FPT2_RPkt(btree,y,D);
perc=100;
if (1)
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
end

img = IPT2_RPkt(btree,tcoef,D);


figure;imagesc(y);title('y');
figure;imagesc(img);title('img')
figure;imagesc(y-img);title('y-img');

return;

close all;
i=10,

for j=1:10:80
  z=zeros(size(tcoef));
  z(i,j)=tcoef(i,j);
  figure;
  wigb(IPT2_RPkt(btree,z,D),1);
end