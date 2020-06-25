% Experiment with pef
close all
clear

plotf=0;
lp=0;

%load /home/dtrad/work3/hyperb2.bin
%data=reshape(hyperb2,512,64);

load /home/dtrad/work3/cdp1000.bin
data=reshape(cdp1000,1751,183);
data=data(513:1024,93:183);

[n m]=size(data);

dataf=fft(data);
figure(1); wigb(data,2);title('data')
figure(2); wigb(fft(dataf.').');title('fk spectrum');


datapf=dataf(1,:);datapf=datapf(:).';
for f=2:n/2-lp;
x=(dataf(f,:));
figure(3)


th=ar(x,10);
%th=armax(x(:),[2 1]);
u=rand(m,1);
y1=predict([x(:) 0*u(:)],th,1);
if (plotf==1)
  subplot(311);plot(real(x));title('x');
  subplot(312);plot(real(y1));title('y1');
  subplot(313);plot(real(y1(:)-x(:)));title('res');figure(gcf)
end
datapf=[datapf;y1(:).'];

end
datapf=[datapf;zeros(lp,m)];
datap=ifft(duplic(datapf));
datap(:,1)=0.001;
figure(4)
wigb(datap,2);title('predicted data');
figure(5)
wigb(data-datap,2);title('residuals');

