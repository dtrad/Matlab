function [data3]=fxprediction2
%example fxprediction interpolation for regular aliased data;
noise = 0.1;
over = 2;
load lineregnoalias
[nt nh]=size(data);

% Insert zero traces in original data;
zerotrace = zeros(nt,1);
data2=data(:,1);
for i=1:nh-1
    data2(:,2*i)=zerotrace;
    data2(:,2*i+1)=data(:,i+1);
end
data2(:,2*nh)=zerotrace;

% generate h2 regular from h 
dh = (h(end)-h(1))/(length(h)-1);
h2=h(1):dh/over:h(end);

%datan=data+(rand(size(data))-0.5)*1;
data3=filter(data2,t,h);

figure,wigb(data2);title('orig')
figure,wigb(data3);title('clean');
% Plot spectra
figure, imagesc(abs(fftshift(fft2(data))));title('spectrum for data');
figure; imagesc(abs(fftshift(fft2(data2))));title('spectrum for zero int');
figure; imagesc(abs(fftshift(fft2(data))));title('spectrum for Spitz');
return;

function [model]=filter(d,t,h)
[nt nh]=size(d);
D=fft(d,2*nt);
DF=zeros(size(D));
dt=t(2)-t(1);
fs=1/dt;nf=nt;
w=2*pi*(0:nf-1)*fs/nt;
DH=D(1:nf,:);
MH=zeros(nf,nh);

for j=1:nf
    MH(j,:) = filterfreq(DH(j,:));
end
M=duplic(MH);

model = ifft(M);model=model(1:nt,:);

return

function [m]=filterfreq(d)
% create pred filter for known data
nh=length(d);nh3=nh+3;
d1=[d(1:2:nh) 0 0];
d2=[0 d(1:2:nh) 0];
d3=[0 0 d(1:2:nh)];
A=[d1(:) d2(:) d3(:)];
[nd np]=size(A);
d4=[d(3:2:nh) 0 0 0];
p=(inv(A'*A+eye(np)*1e-7))*A'*d4(:);

d1=[d(1:nh) 0 0];
d2=[0 d(1:nh) 0];
d3=[0 0 d(1:nh)];
A=[d1(:) d2(:) d3(:)];

% apply for all data
d5=A*(p.^0.5);
m=[d(1) d5(:).'];
m=m(1:nh);
return;

% Plot data
figure, wigb(data);title('Original data');
figure; wigb(data2);title('zero int');
figure; wigb(data);title('Spitz');


