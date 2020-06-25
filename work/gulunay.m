%example Gulunay interpolation for regular aliased data;
noise = 20;

load linereg
[nt nh]=size(data);

% create padded data for operator numerator;
data2=[data;zeros(size(data))];
data3=[data2,zeros(size(data2))];
clear data2;

% kill every second for operator denominator;
data4=data3;
[nt4 nh4]=size(data4);
for i=1:2:nh4
    data4(:,i)=zeros(nt4,1);
end

% Create operator

opn=abs(fft2(data3)); % numerator;
opd=abs(fft2(data4)); % denominator;
op = opn./(opd+noise); % kill operator;

% need to take only first half of frequencies;
op2 = fftshift(op);
nn = (nt+1)/2;
op3 = ifftshift(op2(nn:end-nn,1:end));

clear opn;
clear opd;
clear op;
clear op2;

% Insert zero traces in original data;
zerotrace = zeros(nt,1);
data5=data(:,1);
for i=1:nh-1
    data5(:,2*i)=zerotrace;
    data5(:,2*i+1)=data(:,i+1);
end
data5(:,2*nh)=zerotrace;

% multiply data spectrum for operator;

data6=ifft2(fft2(data5).*op3);

% Plot data
figure, wigb(data);title('Original data');
figure; wigb(data5);title('zero int');
figure; wigb(data6);title('Gulunay');

% Plot spectra
figure, imagesc(abs(fftshift(fft2(data))));title('spectrum for data');
figure; imagesc(abs(fftshift(fft2(data5))));title('spectrum for zero int');
figure; imagesc(abs(fftshift(fft2(data6))));title('spectrum for gulunay');

% Plot intermediate operators;
figure,wigb(data3);title('data for operator numerator');
figure,wigb(data4);title ('data for operator denominator');
