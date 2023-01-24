function mwni_test2(doplot,method,firstfig)

% SSA interpolation for 2d
% irregular data testing 
ninterp=1;

if (nargin < 3) firstfig = 1; end;
if (nargin < 2) method = 'reg';end
if (nargin < 1) doplot = 1;method='reg';end
if (firstfig==1) close all;end

method=method(1:3);

if (method == 'reg') 
    load data2.mat;
else
    load datairreg2.mat
end

nz=10;% zero padding in Kx

data=data(1:512,:);
t=t(1:512);
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

dataorig=data;
if (1==1)
    dataorig=[data(:,1:20) zeros(nt,4) data(:,25:end)];
    %Introduce a gap.
    hh=[h(1:20) h(25:end)];
    h=hh;
    data=[dataorig(:,1:20) dataorig(:,25:end)];
    %%%%%%%%%%%%%%%%%%%%%%
end

% create output axis
nh=length(h)
dx = dh/ninterp;
x=h(1):dx:h(end);
nx=length(x);
%%%%%%%%%%%%%%%%%%%%%%


nk=nx+nz;
map=mapping(h,x);
model=upsample(data,map,nx);
size(model)

    
M=fft(model,nt); 
fs=1/dt;nf=nt/2;
df=1./(dt*nt);
w=2*pi*(0:nf-1)*df;
maxfreqHz=80;
maxfreq=maxfreqHz/df;
MH=M(1:nf,:);     % data in x,f

weights=ones(1,nk);

NITER=zeros(size(1:maxfreq));
[weights]=band_limit(weights,nk);
for f=1:maxfreq
    MMH=embedding(MH(f,:),10);
    MMHf=svdfiltering(MMH,5);
    MH(f,:)=stackantidiagonal(MMHf).';
end
MH(maxfreq+1:end,:)=0;
M=duplic(MH);
model=ifft(M);model=real(model(1:nt,:));
size(model)
index1=1:length(hi);
index2=1:1/ninterp:length(hi);


xi=interp1(index1,hi,index2,'spline');


ifig=firstfig;
if (doplot>1)
    figure(ifig+0),wigb(data,1,hi,t);
    if(doplot >2) print -dpng mwni_orig.png;end

    figure(ifig+1),wigb(model,ninterp,xi,t);
    if(doplot >2) print -dpng mwni_int.png;end

    figure(ifig+2),imagesc(fftshift(abs(fft2(data))));
    if(doplot >2) print -dpng mwni_origspec.png;end

    figure(ifig+3),imagesc(fftshift(abs(fft2(model))));
    if(doplot >2) print -dpng mwni_intspec.png;end

    
    datap=downsample(model,map,nh);
    factor=(max(max(data-datap)))/(max(max(data)));
    figure(ifig+4),wigb(data-datap,factor,hi,t);
    if(doplot >2) print -dpng mwni_res.png;end

end
figure;imagesc(dataorig);title('DATA')
figure;imagesc(model);title('FILTERED DATA');
save model.mat
return

function [M]=embedding(D,k)
for i=1:k
    M(:,i)=wshift(1,D(:),i-1);
end
return

function [Mf]=svdfiltering(M,p)
[U,S,V]=svd(M);
SS=zeros(size(S));    
for i=1:p 
    SS(i,i)=S(i,i);
end
Mf=U*SS*V';
return;
function [out]=stackantidiagonal(a)
%out=a(:,1);
[m,n] = size(a);
idx = hankel(1:m,m:(n-1)+m);
out = accumarray(idx(:),a(:));
out=out(1:m);
return
function [d]=sampling(M,nh,nk,nx,map,weights)
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



function [out]=duplic(in)
% Given half of the FFT series duplic.m  duplicates the data to recover 
% the full spectrum such that the original time series is real.
% the input can be a single vector or a matrix, with the frequency dimension
% longer than the spatial dimension
%
% [out]=duplic(in)
%
%  Daniel Trad, UBC- 30/03/98

[NF NP]=size(in);if NF<NP in=in.';end
[NF NP]=size(in);
 m=1:NF;
 j=1:NP;
 out=zeros(2*NF,NP);
 out(m,j)=in(m,j);

% Nyquist Freq is not in the data so that it must be set 
% either = 0 or = previous value.
out(NF+1,j)=real(in(NF,j));  % Nyquist

 %out(NF+1,j)=real(out(NF,j));
 %out(NF+1,j)=(out(NF,j));   
 
 m=1:NF-1;
 out(NF+1+m,j)=conj(in(NF+1-m,j));
 
 %size(out)
 %size(in)
return;
