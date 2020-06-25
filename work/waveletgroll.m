function [yd,y]=waveletgroll(perc)
% Denoise by thresholding in the ridgelet domain
% testridgeden(perc,noise)
clip=[-0.2,0.2];

close all
t0=129;N=128;
if (nargin<1) perc=0.005;end
% load data
load ozdata0win

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
[yn,wy,twy]=applywt2(y,N,N,perc);
%save wy.mat wy

%load wy.mat
figure(2);imagesc(wy);title('ridgelet coefficients');colorbar;
figure(gcf); 

% Thresholding 
%q=threshold(wy,perc);
%twy=wy;twy(find(abs(wy)<q))=0;

%ttwy=zeros(size(twy));
%ttwy(1:64,1:64)=twy(1:64,1:64);


figure(3);imagesc(twy);title('thresholded ridgelet coefficients');colorbar;
figure(gcf)
% Inverse ridgelet transform
%yn=Inv_FastOrthoRidgeletTrans(twy);
%yn=real(yn);

figure(4);imagesc(yn,clip);title('Noise');
figure(gcf)

yd=y-yn;
figure(5);imagesc(yd,clip);title('Denoised data');
figure(gcf)

return

t=1:N;
figure(6);wigb(y,1,t,t,0.22);title('data');
figure(7);wigb(yd,1,t,t,0.22);title('denoised data');





figure(6);
mx=max(max(abs(y)))*1.2;
axisv=[1 N -mx +mx]


for i=40:2:60
    %subplot(311); plot(x(:,i));axis(axisv);
    subplot(312); plot(y(:,i));axis(axisv);
    subplot(313); plot(x2(:,i));axis(axisv);
    pause
end
figure(2);
figure(5);
