function mwni_test(doplot)
% mwni interpolation for 2d
if (nargin == 0) doplot = 0;end
close all;
load data4.mat
nxz=64;% zero padding in Kx
nyz=32;

data=data(1:550,:,:);
data=data+0*rand(size(data));
t=t(1:550);
nt=length(t);dt=t(2)-t(1);
dhx = hx(2)-hx(1);
dhy = hy(2)-hy(1);

%Introduce a gap.
if (0==1)
    hhx=[hx(1:20) hx(25:end)];
    hhy=[hy(1:5) hy(6:end)];

    hx=hhx;
    hy=hhy;
end

data2=[data(:,1:20,:) data(:,25:end,:)];
data=data2;clear data2;
%%%%%%%%%%%%%%%%%%%%%%

% create output axis
nhx=length(hx)
dx = dhx/2;
x=hx(1):dx:hx(end);
nx=length(x);

nhy=length(hy)
dy = dhy/2;
y=hy(1):dy:hy(end);
ny=length(y);


%%%%%%%%%%%%%%%%%%%%%%


nkx=nx+nxz;
nky=ny+nyz;

mapx=mapping(hx,x);
mapy=mapping(hy,y);

model=resample(data,mapx,nx,mapy,ny);
size(model)
%return;
    
D=fft(data,2*nt);
fs=1/dt;nf=nt;
w=2*pi*(0:nf-1)*fs/nt;
DH=D(1:nf,:);     % data in x,f
MH=zeros(nf,nx); % model in x,f
weights=ones(1,nk);
test=dotproduct_test([DH(1,:) zeros(1,nk)],fft(MH(1,:),nk),weights,nh,nk,nx,map);
display('test = ');test


for f=2:nf-2
    [weights]=band_limit(weights,nk);
    [MH(f,:),weights]=wtcgls_mwni(DH(f,:),nh,nk,nx,map,weights,100);
    if (doplot==1)
        subplot(311);plot(abs(fft(DH(f,:))));
        subplot(312);plot(abs(fft(MH(f,:))));
        subplot(313);plot(weights);figure(gcf);pause;
    end
end
M=duplic(MH);
model=ifft(M);model=model(1:nt,:);
size(model)

if (doplot==2)
    figure(1),wigb(data,1,h,t);
    figure(2),wigb(model,2,x,t);

    figure(3),imagesc(fftshift(abs(fft2(data))));
    figure(4),imagesc(fftshift(abs(fft2(model))));
end

if (doplot==3)
    figure(1),wigb(data,1,h,t);
    print -dpng mwni_orig.png
    
    figure(2),wigb(model,2,x,t);
    print -dpng mwni_int.png

    figure(3),imagesc(fftshift(abs(fft2(data))));
    print -dpng mwni_origspec.png
    
    figure(4),imagesc(fftshift(abs(fft2(model))));
    print -dpng mwni_intspec.png
    
end


return


function [d]=sampling(M,nh,nk,nx,map,weights);
d=zeros(nhx,nhy,1);
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
    
function [m]=resample(d,mapx,nx,mapy,ny)
[nt,nhx,nhy]=size(d);
m=zeros(nt,nx,ny);
for i=1:nhx
    for j=1:nhy
        if ((mapx(i)>0)&&(mapy(j)>0))
            m(:,mapx(i),mapy(j))=d(:,i,j);
        end
    end
end
return;

function [map]=mapping(h,x)
nx=length(x);
nh=length(h);
map=zeros(nh,1);
for i=1:nh
    for j=1:nx
        if (round(h(i))==round(x(j))) map(i)=j;end
    end
end
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
    rms=sqrt(e/length(e));
    if (rms < 1) break;end
end
X=X(:).*w(:);
x=ifft(X);x=x(1:nx);

%update weights,
XX=fft(x,nk);

w=abs(XX);w=smoothing(w,5);
%w=ones(size(w));

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


