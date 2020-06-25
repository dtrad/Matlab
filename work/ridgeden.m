% Denoise by thresholding in the ridgelet domain
close all
clear all

% Generate data
N=128;
x0=zeros(N,N);
x0(50,50)=1;x0(70,70)=1;s
x=(real(FastSlantStack(x0)));
x=x((N/2)+1:(N/2)+N,(N/2)+1:(N/2)+N);
w=ricker(80,0.004);

plot(w);
length(w)
% convolve
wz=[w(:);zeros(N-length(w),1)];
w2=ones(N,1)*wz(:).';
w2=w2';
xx=real(ifft(fft(w2).*fft(x)));

x=xx./norm(xx(:));
figure(1);imagesc(x);title('clean data');colorbar;
figure(gcf)


y=x+0.05*randn(N);
figure(2);imagesc(y);title('noisy data');colorbar;
figure(gcf)

% Apply ridgelet transform
wy=FastOrthoRidgeletTransform(y);
figure(3);imagesc(wy);title('ridgelet coefficients');colorbar;
figure(gcf);

% Thresholding 
% Define the threshold such that only 100 coefficients remain
total=2*N*2*N;
kept=round((total)*0.005)
message=sprintf('number of coefficients to keep=%d from %d or %f percent',...
    kept,total,kept/total*100)



ss=sort((abs(wy(:))));q=ss(length(ss)-kept);
twy=wy;twy(find(abs(wy)<q))=0;
figure(4);imagesc(twy);title('thresholded ridgelet coefficients');colorbar;
figure(gcf)
% Inverse ridgelet transform
x2=Inv_FastOrthoRidgeletTrans(twy);
x2=real(x2);
figure(5);imagesc(x2);title('Denoised data');colorbar;
figure(gcf)

figure;
mx=max(max(abs(y)))*1.2;

axisv=[1 N -mx +mx]
for i=10:30
    subplot(311); plot(x(:,i));axis(axisv);
    subplot(312); plot(y(:,i));axis(axisv);
    subplot(313); plot(x2(:,i));axis(axisv);
    pause
end
figure(2);
figure(5);
