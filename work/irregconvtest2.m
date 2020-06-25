function irregconvtest(method)

xmin = -20;
xmax = 21;
interval = .1;
xtest=[xmin:interval:xmax];
ntest = length(xtest);
q=20;
f=2;
variation = 0.2;
x=xtest;
pert = rand(size(x))-0.5;
x=x+variation*pert;
freq1=0.1;
freq2=0.2;
y=sin(2*pi*x*freq1)+sin(2*pi*x*freq2);
ytest=sin(2*pi*xtest*freq1)+sin(2*pi*xtest*freq2);

if (method == 'std')
    [xr,zr]=irregconv(x,y,q);
elseif (method == 'dui')
    [xr,zr]=duijn(x,y,q,f);
end

ztest = conv(ytest,func(xtest));


figure(1)
subplot(221),plot(x,func(x));title('Gaussian');
subplot(222),plot(x,y), title('irregularly sampled signal');
subplot(223), plot(1:length(zr),zr);title('irreg convolution g*f')
subplot(224),plot(xtest,ztest(ntest/2+1:ntest+ntest/2));title('regular convolution g*f ');

figure(gcf)
keyboard;
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

function [xr,zr]=irregconv(x,y,q)
% performs discrete conv for a irregular function and an analytical function 
% z(t)=sum(f(tau)*g(t-tau))

%indexes for regular axis
nx=length(x);
xr = [min(x):1:max(x)];
nxr = length(xr);
deltax = xr(2)-xr(1);

% padding q points for truncation
minx = -floor((q+1)/2);
maxx = floor((q+1)/2) + nxr;
xc = minx:1:maxx;
xc = xc*deltax;

nc = length(xc);
zr = zeros(1,nc);

% save deltas and func for visualization
deltaMatrix = zeros(nc,nx);
funcMatrix = zeros(nc,nx);

for i=1:nc
    for j=1:nx
        deltaMatrix(i,j) = xc(i)-x(j);
        delta = deltaMatrix(i,j);
        ratio = abs(delta/deltax); % truncate filter
        if (ratio < (q/2)) funcMatrix(i,j)= func(delta);
        else funcMatrix(i,j) = 0;
        end
        % which is the appropriate dx?
        dx=xc(i) - x(j);        
        zr(i)=zr(i)+dx*y(j)*funcMatrix(i,j);
    end
end

zr = zr(1:end-q-1);
size(zr)

return

function [xc,p2]=duijn(x,y,q,f)
% performs discrete conv for a irregular function and an analytical function 
% z(t)=sum(f(tau)*g(t-tau))

nx=length(x);
%normalize axis to [0-1]
xn = x - x(1);
xn = xn/xn(nx);

Nc = f*nx; % oversampling
deltax = 1/Nc;
xc = -floor((q+1)/2)+1:Nc+floor((q+1)/2)-1;
nc = length(xc);
p2 = zeros(1,nc+floor(q/2));
%keyboard
funcMatrix=zeros(nc,nx);
deltaMinv = zeros(1,nx);
deltaMaxv = zeros(1,nx);
for j=1:nx
    ndelta = floor(xn(j)/deltax);
    deltaMin = ndelta-q/2+1;
    deltaMax = ndelta+q/2;
    deltaMinv(j)=deltaMin;
    deltaMaxv(j)=deltaMax;
    for k=deltaMin:deltaMax
        delta = k*deltax - xn(j);
        index = k + floor(q/2);
        signinc = sign(delta);
        %signinc = 1;
        p2(index)=p2(index)+signinc*deltax*y(j)*func(delta);
    end
end
p2=p2(1:end-q/2);
size(p2)
size(xc)

%plot(p2),figure(gcf);
%keyboard

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
    y(i)=exp(-x(i)*x(i)/10);
end
return

function [xn]=normaxis(x);
nx=length(x);
xn = x - x(1);
xn = xn/xn(nx);
return;