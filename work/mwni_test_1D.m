function [data,h]=mwni_test_1D(doplot,method,firstfig)

% mwni interpolation for 1d
% irregular data testing 
ninterp=3;

if (nargin < 3) firstfig = 1; end;
if (nargin < 2) method = 'reg';end
if (nargin < 1) doplot = 0;method='reg';end
if (firstfig==1) close all;end

method=method(1:3);
A1=1;
A2=2;
f1 = 0.05;
f2 = 0.2;
nh=64;
h = 0: nh-1;
if (method == 'reg') 
    display('regular');
else
    h=h+1*(rand(1,nh)-0.5);
end
data = A1*sin(2*pi*f1*h)+A2*sin(2*pi*f2*h);

nz=nh*ninterp;% zero padding in Kx
nz=10;
data=data(:).';


if (0==1) 
    %Introduce a gap.
    hh=[h(1:20) h(25:end)];
    h=hh;
    data2=[data(1:20) data(25:end)];
    data=data2;
    %%%%%%%%%%%%%%%%%%%%%%
end

figure(1);plot(h,data,'-o');title('input data')
figure(2);plot(fftshift(abs(fft(data))));

% generate h regular from h irregular
% (even if h is already regular
dh = (h(end)-h(1))/(length(h)-1);
hi=h;
h=hi(1):dh:hi(end);
if (length(h) ~= length(hi))
    display('error');
    return;
end


% create output axis
nh=length(h)
dx = dh/ninterp;
x=h(1):dx:h(end);
nx=length(x);
%%%%%%%%%%%%%%%%%%%%%%

nk=nx+nz;

% mapping from binned input to regular output
[map,binerror] = mapping(hi,x);
hiorig=hi;
[map,hi,data,nh]=eliminate_rep(map,hi,data);
model=upsample(data,map,nx);
nh
ii=1:nh-1;find((map(ii+1)-map(ii))==0);
size(model)

weights=ones(1,nk);
%test=dotproduct_test([DH(1,:) zeros(1,nk)],fft(MH(1,:),nk),weights,nh,nk,nx,map);
%display('test = ');test

[weights]=band_limit(weights,nk);
for ii=1:3;
[model,weights]=wtcgls_mwni(data,nh,nk,nx,map ,weights,50);
end

if (doplot==1)
    figure(1);
    subplot(311);plot(abs(fft(data)));
    subplot(312);plot(abs(fft(model)));
    subplot(313);plot(weights);figure(gcf);pause;
end

display('output size is ');size(model);

% create an irregular axis with same shape as hi
index1=1:length(hiorig);
index2=1:1/ninterp:length(hiorig);
xi=interp1(index1,hiorig,index2,'spline');



%figure(1);plot(hi,real(data),'-o')
%figure(2);plot(x,real(model),'-o')
%figure(3);plot(xi,real(model),'-o')
datafulli = A1*sin(2*pi*f1*xi)+A2*sin(2*pi*f2*xi);
datafullr = A1*sin(2*pi*f1*x)+A2*sin(2*pi*f2*x);

figure(3);plot(xi,datafulli,x,real(model),'o')
figure(4);plot(xi,datafulli,xi,real(model),'o')
figure(5);plot(hi,data,'-o',hi,real(model(map)),'+');title('pred orig');

model2= model;
model2(map)=0;
figure(6);plot(hi,data,'-o',x,real(model2),'+');title('pred new');

%figure(4);plot(xi,datafulli,hi,data,'o')
%keyboard;

ifig=firstfig;
if (doplot>1)
    figure(ifig+0),wigb(data,1,hi,t);
    if(doplot >2) print -dpng mwni_orig.png;end

    figure(ifig+1),wigb(model,ninterp,xi,t);
    if(doplot >2) print -dpng mwni_int.png;end

    figure(ifig+2),wigb(model,ninterp,1:length(xi),t);
    if(doplot >2) print -dpng mwni_intreg.png;end

    figure(ifig+3),imagesc(fftshift(abs(fft2(data))));
    if(doplot >2) print -dpng mwni_origspec.png;end

    figure(ifig+4),imagesc(fftshift(abs(fft2(model))));
    if(doplot >2) print -dpng mwni_intspec.png;end

    
    datap=downsample(model,map,nh);
    factor=(max(max(data-datap)))/(max(max(data)));
    figure(ifig+5),wigb(data-datap,factor,hi,t);
    plot(binerror*1/dh);
    if(doplot >2) print -dpng mwni_res.png;end

end

return

function [d]=sampling(M,nh,nk,nx,map,weights);
d=zeros(nh,1);
weights=weights(:);
M=M(:);
M=M.*weights;
m=ifft(M);m=m(1:nx);
checkshape(d,m,weights);
for i=1:nh
    d(i)=m(map(i));
end
d=d(:);

return;

function [M]=interpolate(d,nh,nk,nx,map,weights)
weights=weights(:);
m=zeros(nx,1);
d=d(:);
weights=weights(:);
checkshape(d,m,weights);
for i=1:nh
    m(map(i))=d(i);
end
M=fft(m,nk)/nk;
M=M(:).*weights(:);


return;

function checkshape(d,x,w)
[n m]=size(d);
if (m>1) error('d is not vector');size(d),end
[n m]=size(x);
if (m>1) error('x is not vector');end
[n m]=size(w);
if (m>1) error('w is not vector');end
return    
    
function [m]=upsample(d,map,nx)
[nt,nh]=size(d);
m=zeros(nt,nx);
for i=1:nh
    if (map(i)>0)
        m(:,map(i))=d(:,i);
    end
end
return;

function [d]=downsample(m,map,nh)
[nt,nx]=size(m);
d=zeros(nt,nh);
for i=1:nh
    if (map(i)>0)
        d(:,i)=m(:,map(i));
    end
end
return;

function [map,binerror]=mapping(h,x)
nx=length(x);
nh=length(h);
dx = x(2)-x(1);
map=zeros(nh,1);
for i=1:nh
        map(i)=round((h(i) - x(1))/dx) + 1;
        binerror(i)=h(i)-x(1)+map(i)*dx;
end
return

function [map,binerror]=binning(h,x)
nx=length(x);
nh=length(h);
map=zeros(nh,1);
binerror=zeros(size(map));
for i=1:nh
    for j=1:nx-1
        if ((h(i)>=x(j))&(h(i)<x(j+1))) 
            map(i)=j;
            binerror(i)=h(i)-x(j);
            %need extra check because of possible numerical error
            if (binerror(i) > 0) 
                binerrornext=h(i)-x(j+1);
                if ((abs(binerrornext))< (abs(binerror(i))))
                    map(i)=j+1;
                    binerror(i)=binerrornext;
                end
            end
            break;
        end
    end
end
if (h(end) == x(end)) map(nh) = nx;end
return


function [x,w]=wtcgls_mwni(d,nh,nk,nx,map,w,niter)
eps=1e-3;
x=zeros(nx,1);
X=zeros(nk,1);
d=d(:);
y=[d;X];

y0=sampling(X,nh,nk,nx,map,w);
y0(nh+1:nh+nk)=eps*X(1:nk);
rr=y-y0;

energy=rr'*rr;
e_old=energy;

G=interpolate(rr,nh,nk,nx,map,w);
G=G+eps*rr(nh+1:nh+nk);

S=G;
gammam=G'*G;
rms=0;

for iter=1:niter
    ss=sampling(S,nh,nk,nx,map,w);
    ss(nh+1:nh+nk)=eps*S(1:nk);
    den=ss'*ss;
    %denp(iter)=den;plot(denp);figure(gcf);pause
    %if (abs(den)<1e-7) keyboard;end
    alpha=gammam/den;
    X=X+alpha*S;
    rr=rr-alpha*ss;
    G=interpolate(rr,nh,nk,nx,map,w);
    G=G+eps*rr(nh+1:nh+nk);
    gamma=G'*G;
    %if (abs(gammam)<eps) return;end
    beta=gamma/gammam;
    gammam=gamma;
    S=G+beta*S;
    e=rr'*rr;
    e_old=e;
    rmsold =rms;
    rms=sqrt(e/length(e));
    if (rms == rmsold) break;end
    rmsv(iter)=rms;
    if (rms < .1) break;end
end
iter
X=X(:).*w(:);
x=ifft(X);x=x(1:nx);

%update weights,
XX=fft(x,nk);

w=abs(XX);w=smoothing(w,5);
%w=ones(size(w));

figure;plot(rmsv);

return

function [w]=band_limit(w,nk)
li=round(nk/4+1);
ui=round(nk*3/4);
w(li:ui)=0;
return;

function [y]=smoothing(x,nl)
y=x;nx=length(x);
nl2=(nl-1)/2;
sum=0;
for i=1:nl
    sum = sum + x(i);
end
for i=nl2+2:nx-nl2
    sum = sum - x(i-nl2-1) + x(i+nl2);
    y(i)=sum/nl;
end
return;

function [test]=dotproduct_test(d,m,w,nh,nk,nx,map);
D1=rand(size(d))+i*rand(size(d));
M1=rand(size(m))+i*rand(size(m));
D1=D1(:);
M1=M1(:);
weights=rand(size(w));
[D2]=sampling(M1,nh,nk,nx,map,weights);
D2(nh+1:nh+nk)=eps*M1(1:nk);

[M2]=interpolate(D1,nh,nk,nx,map,weights);
M2=M2+eps*D1(nh+1:nh+nk);

dot1=D1(:)'*D2(:);
dot2=M1(:)'*M2(:);

 test=abs(dot1/dot2);
 return;

 function [map2,h2,data2,nh2]=eliminate_rep(map,h,data);
 nh=length(map); 
 map2=map(1);
 h2=h(1);
 data2=data(1);
 for i=2:nh
     if (map(i) ~= map(i-1)) 
         map2=[map2;map(i)];
         h2=[h2;h(i)];
         data2=[data2;data(i)];
     end
 end
 nh2=length(map2);
 display('removed ');
 nh - nh2
 return;
