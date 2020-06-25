
function filter_yuval

t=1:512;x=(1*sin(2*pi*0.1*t)+0.8*sin(2*pi*0.2*t));
subplot(221)
plot(x);

subplot(222);
plot(abs(fft(x)))
y=real(filter_yuval_enemies(x));

subplot(223);
plot(y);

subplot(224);
plot(abs(fft(y)));
figure(gcf)
return

function [y]=filter_yuval_enemies(x)
nx=length(x);

xx = [x zeros(size(x))];
nxh = length(xx)/2;

XX=fft(xx);
XXH = XX(1:nxh);
threshold = max(abs(XXH))*0.9

I=find(abs(XXH)>threshold);
XXH(I)=0;
yy=ifft(duplic(XXH));

y=yy(1:nx);

return
