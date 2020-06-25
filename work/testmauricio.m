

clear; close all;


% Regular full sampling;
ny = 16;
nx = 16;


s = ones(ny,nx); 
S = fftshift(abs(fft2(s))); 


figure(1); 

subplot(221); plot_grid(s);  axis tight; 
xlabel('x'); 
ylabel('y'); 
subplot(222); imagesc(S);colormap(1-gray);
xlabel('kx'); 
ylabel('ky'); 


% decimation in x

I=1:2:nx;
s(:,I)=0;
S = fftshift(abs(fft2(s))); 


subplot(223); plot_grid(s);  axis tight; 
xlabel('x'); 
ylabel('y'); 
subplot(224); imagesc(S);colormap(1-gray);
xlabel('kx'); 
ylabel('ky'); 



% decimation in x and y

I=1:2:nx;
J=1:2:ny;
s(:,I)=0;
s(J,:)=0;
S = fftshift(abs(fft2(s))); 


figure(2); 

subplot(221); plot_grid(s);  axis tight; 
xlabel('x'); 
ylabel('y'); 
subplot(222); imagesc(S);colormap(1-gray);
xlabel('kx'); 
ylabel('ky'); 




% random decimation in x and y

s = round(full(sprand(ny,nx,0.5)));
S = fftshift(abs(fft2(s))); 


figure(2); 

subplot(223); plot_grid(s);  axis tight; 
xlabel('x'); 
ylabel('y'); 
subplot(224); imagesc(S);colormap(1-gray);
xlabel('kx'); 
ylabel('ky'); 












