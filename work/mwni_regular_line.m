function [log]=mwni_test2(doplot,method,firstfig)
% 
% Global interpolator: mwni interpolation for 2d
% irregular data testing
% Output is irregular (preserves original data).
% it creates a regular plot by using oversampling / binning
ninterp=4; % set to 6 for oversampling (regularization), 2 for interp.

% Parameters: 
%           doplot = 1: plot iteration by iteration
%           doplot = 2: normal plots on screen
%           doplot = 3: normal plots on screen and disk
%           method = 'reg', 'ireg'
%           firstfig = 1; use to preserve previous plots on screen.
% Daniel Trad; Veritas, March, 2006 (last May 4, 06)
%                             
if (nargin < 3) firstfig = 1; end;
if (nargin < 2) method = 'reg';end
if (nargin < 1) doplot = 0;method='reg';end
if (firstfig==1) close all;end

method=method(1:3);

if (method == 'reg') 
    load data2.mat;
else
    load lineirreg.mat
end

nz=5;% zero padding in Kx

data=data(1:550,:);
t=t(1:550);
nt=length(t);dt=t(2)-t(1);

% generate h regular from h irregular
% (even if h is already regular
dh = (h(end)-h(1))/(length(h)-1);
hi=h;
h=hi(1):dh:hi(end);
if (length(h) ~= length(hi))
    display('error');
    return;
end


if (1==0) 
    %Introduce a gap.
    hh=[h(1:20) h(25:end)];
    h=hh;
    data2=[data(:,1:20) data(:,25:end)];
    data=data2;
    %%%%%%%%%%%%%%%%%%%%%%
end

% create output axis
nh=length(h)
dx = dh/ninterp;
x=h(1):dx:h(end);
nx=length(x);
%%%%%%%%%%%%%%%%%%%%%%
%irregular axis for output
index1=1:length(hi);
index2=1:1/ninterp:length(hi);
xi=interp1(index1,hi,index2,'spline');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nk=nx+nz;
map=mapping(h,x);
model=upsample(data,map,nx);
size(model)

D=fft(data,2*nt);
fs=1/dt;nf=nt;
w=2*pi*(0:nf-1)*fs/nt;
DH=D(1:nf,:);     % data in x,f
MH=zeros(nf,nx); % model in x,f
weights=ones(1,nk);
test=dotproduct_test([DH(1,:) zeros(1,nk)],fft(MH(1,:),nk),weights,nh,nk,nx,map);
display('test = ');test
niter=5;
[weights]=band_limit(weights,nk,ninterp);
log=zeros(nf,2);
for f=5:nf-30   
    [MH(f,:),weights,iter]=wtcgls_mwni(DH(f,:),nh,nk,nx,map,weights,niter);    
    if (doplot==1)
        subplot(311);plot(abs(fft(DH(f,:))));
        subplot(312);plot(abs(fft(MH(f,:))));
        subplot(313);plot(weights);figure(gcf);pause;
    end
    log(f,:)=[f iter];
end
log=log(5:nf-30,:);

M=duplic(MH);
model=ifft(M);model=model(1:nt,:);
size(model)
ifig=firstfig;
if (doplot>1)
    figure(ifig+0),wigb(data,1,hi,t);
    if(doplot >2) print -dpng mwni_orig.png;end

    figure(ifig+1),wigb(model,ninterp,xi,t);
    if(doplot >2) print -dpng mwni_int.png;end
    
    %define binned axis with space 25;
    [xb,index]=binning(xi,25);
    figure(ifig+2),wigb(model(:,index),1,xb,t);
    if(doplot >2) print -dpng mwni_reg.png;end

    figure(ifig+3),imagesc(fftshift(abs(fft2(data))));
    if(doplot >2) print -dpng mwni_origspec.png;end

    figure(ifig+4),imagesc(fftshift(abs(fft2(model))));
    if(doplot >2) print -dpng mwni_intspec.png;end

    
    datap=downsample(model,map,nh);
    factor=(max(max(data-datap)))/(max(max(data)));
    figure(ifig+5),wigb(data-datap,factor,hi,t);
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


function [x,w,iter]=wtcgls_mwni(d,nh,nk,nx,map,w,niter)
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
    if (abs(den)<1e-14) 
        break;
    end
    alpha=gammam/den;
    X=X+alpha*S;
    rr=rr-alpha*ss;
    G=interpolate(rr,nh,nk,nx,map,w);
    G=G+eps*rr(nh+1:nh+nk);
    gamma=G'*G;
    if (abs(gammam)<1e-4) 
        break;
    end
    beta=gamma/gammam;
    gammam=gamma;
    S=G+beta*S;
    e=rr'*rr;
    e_old=e;
    rms=sqrt(e/length(e));
    %if (rms < 1) break;end
end
X=X(:).*w(:);
x=ifft(X);x=x(1:nx);

%update weights,
XX=fft(x,nk);

w=abs(XX);w=smoothing(w,5);
%w=ones(size(w));

return

function [w2]=band_limit(w,nk,ninterp)
if (ninterp==1) w2=w;return;end
li=round(nk/4+2);
ui=round(nk*3/4);
w2=w;
w2(li:ui)=0;
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
 
function [x,index]=binning(xi,dx)
nx=length(xi);

j=1;
index(j)=1;
x(1)=round(xi(1)/dx)*dx;
for i=2:nx
    temp=round(xi(i)/dx)*dx;
    if (temp ~= x(j))
        j=j+1;
        x(j)=temp;
        index(j)=i;
    end
end
return;
        