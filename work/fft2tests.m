clear 
close all

load data
[nt nh]=size(data);
data(:,1:nh/2)=0;

fig=1;

figure(fig);fig=fig+1;
imagesc(abs(fft2(data)));colorbar;title('original')
figure(fig);fig=fig+1;
imagesc(fftshift(abs(fft2(data))));colorbar;title('fftshift original')


%data2=zeros(size(data));
data2=data;

%center the spectrum in w
for i=1:2:nt, 
  data2(i,:)=data(i,:)*(-1);
end
%center the spectrum in kx
for i=1:2:nh, 
  data2(:,i)=data2(:,i)*(-1);
end



figure(fig);fig=fig+1;
imagesc(abs(fft2(data2)));colorbar;title('shifted spectrum');

