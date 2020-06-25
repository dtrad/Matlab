function irregconvtest
x=[-14 -13.5 -13 -12 -10 -9 -8.4 -7 -6 -5.1 -4 -3 -2.8 -1 0 0.5 1 ...
   1.3 2 2.8  4 4.4 5 7 8 9 10 11 12 13.5 14 15 18 21 23];
x=rand(64,1);x=sort(x-0.5);x=x*128;x=round(x);
%x=-64:64;
nx=length(x);
%x'
xtest=[min(x):1:max(x)];
ntest = length(xtest);
%x=xtest;
freq1=0.1;
freq2=0.2;
y=sin(2*pi*x*freq1)+sin(2*pi*x*freq2);
ytest=sin(2*pi*xtest*freq1)+sin(2*pi*xtest*freq2);
[xr,zr]=irregconv(x,y);
ztest = conv(ytest,func(xtest),'same');
zitest = conv(y,func(x),'same');

figure(1)
subplot(221),plot(x,func(x));title('Gaussian');
%subplot(222),plot(x,y), title('irregularly sampled signal');
subplot(222),plot(x,zitest), title('ignored irreg conv signal');
subplot(223), plot(xr,zr);title('irreg convolution g*f')
subplot(224),plot(xtest,ztest);title('regular convolution g*f ');
%subplot(224),plot(xtest,ztest(round(ntest/2+1):round(ntest+ntest/2)));title('regular convolution g*f ');

figure(gcf)

return;

%deconv
yr=deconi(xr,zr);

figure(2)

subplot(222),plot(xr,yr,x,y);title('reg vs irreg');
subplot(223),plot(xr,zr,x,y);title('filt reg vs irreg')
subplot(224),plot(xr,zr,xr,yr);title('conv vs deconv');

figure(gcf);

figure(3)
subplot(221),plot(xr,abs(fft(yr)));title('spectrum deconv')
subplot(222),plot(xr,abs(fft(zr)));title('spectrum filtered');
subplot(223),plot(xr,abs(fft(func(xr))));title('spectrum deconv')
subplot(224),plot(xtest,abs(fft(ytest)));title('spectrum orig')

return

function [xr,z]=irregconv(x,y)
% performs discrete conv for a irregular function and an analytical function 
% z(t)=sum(f(tau)*g(t-tau))

nx=length(x);
xr = [min(x):1:max(x)];
nxr = length(xr);
z=zeros(1,nxr);

for i=1:nxr
    z(i)=0;
    for j=1:nx
        %z(i)=z(i)+(xr(i)-x(j))*y(j)*func(xr(i)-x(j));
        z(i)=z(i)+y(j)*func(xr(i)-x(j));%*(xr(i)-x(j));
    end
end

%z=z./(nxr/nx);

return

function [yr]=deconi(x,y)
f=func(x);
nf =length(f);

%zero padding
nz = nf*2;
F=fft(f,nz);
Y=fft(y,nz);

YR=(conj(F).*Y)./max((conj(F).*F),1e-5);
nh=length(YR)/2
YRH=YR(1:nh);
YRH(nh-20:nh)=0;
YR=duplic(YRH).';
yr=real(ifft(YR));
yr=yr(1:nf);

return

% work with positive freqs only
FH=F(1:nf);
YH=Y(1:nf);
 
YRH=(conj(FH).*YH)./max((conj(FH).*FH),1e-3);
yr=real(ifft(duplic(YRH)));

%truncate zero padding
yr=yr(1:nf);
yr=yr.';



return


function [y]=func(x)
for i=1:length(x)
    y(i)=exp(-x(i)*x(i)/5);
end
return