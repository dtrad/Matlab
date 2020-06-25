function [yd,y]=ridgefilter(y,perc)
% Denoise by thresholding in the ridgelet domain
% testridgeden(perc,noise)
if (nargin<2) perc=0.005;end
[N M]=size(y)
figure(1)
imagesc(y);title('data');

% Apply ridgelet transform
wy=FastOrthoRidgeletTransform(y);
figure(2);imagesc(wy);title('ridgelet coefficients');colorbar;
figure(gcf); 

% Thresholding 
q=threshold(wy,perc);
twy=wy;twy(find(abs(wy)<q))=0;
figure(3);imagesc(twy);title('thresholded ridgelet coefficients');colorbar;
figure(gcf)
% Inverse ridgelet transform
yn=Inv_FastOrthoRidgeletTrans(twy);
yn=real(yn);

figure(4);imagesc(yn);title('Noise');
figure(gcf)

yd=y-yn;
figure(5);imagesc(y-yn);title('Denoised data');
figure(gcf)

t=1:N;
figure(6);wigb(y,1,t,t,0.22);title('data');
figure(7);wigb(yd,1,t,t,0.22);title('denoised data');

return
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
