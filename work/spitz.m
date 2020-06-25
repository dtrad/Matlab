%example Spitz interpolation for regular aliased data;
noise = 20;

load linereg
[nt nh]=size(data);

% Insert zero traces in original data;
zerotrace = zeros(nt,1);
data2=data(:,1);
for i=1:nh-1
    data2(:,2*i)=zerotrace;
    data2(:,2*i+1)=data(:,i+1);
end
data2(:,2*nh)=zerotrace;

data3=filter(data2);

return;

function [model]=filter(d,t,h)
[nt nh]=size(d);
D=fft(d,2*nt);
DF=zeros(size(D));
dt=t(2)-t(1);
fs=1/dt;
w=2*pi*(0:nf-1)*fs/nt;
DH=D(1:nf,:);
MH=zeros(nf,nh);

for j=1:nf
    MH(j,:) = filterfreq(DH(j,:);
end
M=duplic(MH);

model = ifft(M);model=model(1:nt,:);

return

function [m]=filterfreq(d)
m=d;
return;

% Plot data
figure, wigb(data);title('Original data');
figure; wigb(data2);title('zero int');
figure; wigb(data);title('Spitz');

% Plot spectra
figure, imagesc(abs(fftshift(fft2(data))));title('spectrum for data');
figure; imagesc(abs(fftshift(fft2(data2))));title('spectrum for zero int');
figure; imagesc(abs(fftshift(fft2(data))));title('spectrum for Spitz');
