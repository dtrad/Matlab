function [data,h,F,weights]=mwni_dft_test_1D(doplot,method,dftop,firstfig)

% mwni interpolation for 1d, comparing DFT and FFT
% regular and irregular data testing 
ninterp=3;
printing=1;
ifig=1;%index for incremental figure numbers.

if (nargin < 4) firstfig = 1; end;
if (nargin < 3) dftop  = 'yes';end
if (nargin < 2) method = 'reg';end
if (nargin < 1) doplot = 0;end
if (firstfig==1) close all;end
method = method(1:3);
dftop = dftop(1:2);

if (dftop == 'no') 
    dftop='no ';
    ft='FFT';
elseif (dftop == 'ye')
    dftop='yes';
    ft='DFT';
end

display(sprintf('h axis is %s and DFT is %s ',method,dftop));
   

% generate or read axis h
nh=64;dh=1;
h = 0: nh-1;
if (method ~= 'reg')
    h=h+1*(rand(1,nh)-0.5);
end

save pp.mat h;
%load pp.mat h;

data = my_function(h);
freqdata = freqaxis(dh,nh);
hfine=0:nh*10;hfine=hfine/10;
datafine = my_function(hfine);

nz=10;
data=data(:).';
maxd = max(data);
mind = min(data);
axisv=[0 nh mind maxd];

if (0==1) 
    %Introduce a gap.
    hh=[h(1:20) h(25:end)];
    h=hh;
    data2=[data(1:20) data(25:end)];
    data=data2;
    %%%%%%%%%%%%%%%%%%%%%%
end

figure(ifig);plot(hfine,datafine,h,data,'o');title('input data')
ifig=ifig+1;axis(axisv);

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
nh=length(h);
dx = dh/ninterp;
x=h(1):dx:h(end);
nx=length(x);
freqx=freqaxis(dx,nx);
%%%%%%%%%%%%%%%%%%%%%%



% mapping from binned input to regular output
[map,binerror] = mapping(hi,x);
hiorig=hi; 

%eliminate repeated bins. FFT can handle only one position per bin.
%remove the hi,data values corresponding to repeated bins.
[map,hi,data,nh]=eliminate_rep(map,hi,data);
nhi=length(hi); 
freqdata=freqaxis(dh,nhi);
model=upsample(data,map,nx);
display(sprintf('size of model is %d \n',length(model)));


%test=dotproduct_test([DH(1,:) zeros(1,nk)],fft(MH(1,:),nk),weights,nh,nk,nx,map);
%display('test = ');test

% create an irregular axis with same shape as hi
% this axis is necessary to match the input points exactly,
% which can not be done with a regular axis.
% notice this is a requirement even using the DFT because 
% the least squares algorightm needs to match the irregular input.

% use hiorig in case any of the end points was removed by eliminate_rep
index1=0:length(hiorig)-1;
index2=0:1/ninterp:length(hiorig)-1;
xi=interp1(index1,hiorig,index2,'spline');

if (doplot > 0)
figure(ifig);plot(index1,hiorig,'-o',index2,xi,'-+');ifig=ifig+1;
title('original axis (o) and new axis (+)');
end

nk=nx+nz;
if ( dftop == 'yes' )
    % Distance weights give bad results. Make them one for now.
    wxi=ones(size(1,nx));
    % DFT operator
    % 1) DFT for irregular axis.
    [F,kx]=dftoperator(xi,nz,wxi);
    % 2) DFT for regular axis.
    [FR]=reg_operator(x,kx);
    % Length for the kx axis
    nk=length(kx);    
end
weights=ones(1,nk);
[weights]=band_limit(weights,nk);

nextiter=5;
if (dftop == 'yes')
    for ii=1:nextiter;
        [model,weights,X,iter]=wtcgls_mwni(data,nh,nk,nx,map ,weights,150,F);
    end
else
    for ii=1:nextiter;
        [model,weights,X,iter]=wtcgls_mwni(data,nh,nk,nx,map ,weights,150);
    end
end

if (1)
    % with weigths
    FHF=F*F' + diag(weights);
    [U,S,V]=svd(F'*diag(1./weights+eps));
    figure(ifig);mesh(abs(FHF));ifig=ifig+1;title('F^HF with weights');
    print -dpng FHF_dft_w.png;
    figure(ifig);nnv=length(V(:,1));
    for k=1:1:4;subplot(220+k);plot(1:nnv,real(V(:,k)),1:nnv,imag(V(:,k)));title(sprintf('V(:,%d) vector with weights\n',k));end;;
    print -dpng vectorsV_dft_w.png;
    ifig=ifig+1;
    % without weights
    FHF=F*F';
    [U,S,V]=svd(F');
    figure(ifig);mesh(abs(FHF));ifig=ifig+1;title('F^H F without weights');
    print -dpng FHF_dft_nw.png;
    figure(ifig);
    for k=1:1:4;subplot(220+k);plot(1:nnv,real(V(:,k)),1:nnv,imag(V(:,k)));title(sprintf('V(:,%d) vector no weights\n',k));end;
    print -dpng vectorsV_dft_nw.png;
    ifig=ifig+1;
end

if (dftop == 'yes')
    modelreg = idft(FR,X);
else
    modelreg = model;
end

if (doplot > 0)
    figure(ifig);
    axisvw=[-0.5 0.5 0 max(abs(fft(model)))];
    
    subplot(311);plotspectrum(data,dh);title('fft(data)');axis(axisvw);
    subplot(312);plotspectrum(model,dh/ninterp);title('fft(model)');axis(axisvw);
    subplot(313);plotspectrum(ifft(weights),dh/ninterp);title('weights');axis(axisvw);
    figure(gcf);ifig=ifig+1
end

display(sprintf('output size is %d\n',length(model)));

%figure(1);plot(hi,real(data),'-o')
%figure(2);plot(x,real(model),'-o')
%figure(3);plot(xi,real(model),'-o')
datafulli = my_function(xi);
datafullr = my_function(x);

figure(ifig);plot(x,datafullr,x,real(modelreg),'-o');
title(sprintf('predicted with %s on regular axis',ft));axis(axisv);ifig=ifig+1;
figure(ifig);plot(xi,datafulli,xi,real(model),'-o');
title(sprintf('predicted with %s on irregular axis',ft))
axis(axisv);ifig=ifig+1;
figure(ifig);plot(hi,data,'-o',hi,real(model(map)),'+');title('pred orig');ifig=ifig+1;

if (printing == 1)
if (dftop == 'yes')
    figure(1);print -dpng data_orig_dft.png;
    figure(3);print -dpng spectra_dft.png;
    figure(4);print -dpng data_pred_dft_reg.png;
    figure(5);print -dpng data_pred_dft_irreg.png;
    figure(6);print -dpng data_pred_dft_orig.png;
else 
    figure(1);print -dpng data_orig_fft.png;
    figure(3);print -dpng spectra_fft.png;
    figure(4);print -dpng data_pred_fft_reg.png;
    figure(5);print -dpng data_pred_fft_irreg.png;
    figure(6);print -dpng data_pred_fft_orig.png;
end
end
model2= modelreg;
model2(map)=0;
figure(ifig);plot(hi,data,'-o',x,real(model2),'+');title('pred new');axis(axisv);ifig=ifig+1;

return

function [d]=sampling(M,nh,nk,nx,map,weights,F);
d=zeros(nh,1);
weights=weights(:);
M=M(:);
M=M.*weights;
if ( nargin == 7) % given F operator 
    m = idft(F,M);
else
    m=ifft(M);    % else use FFT
end
    
m=m(1:nx);
checkshape(d,m,weights);
for i=1:nh
    d(i)=m(map(i));
end
d=d(:);

return;

function [M]=interpolate(d,nh,nk,nx,map,weights,F)
weights=weights(:);
m=zeros(nx,1);
d=d(:);
weights=weights(:);
checkshape(d,m,weights);
for i=1:nh
    m(map(i))=d(i);
end
if (nargin == 7)
    M=dft(F,m)/nk;
else
    M=fft(m,nk)/nk;
end
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

function [x,w,X,iter]=wtcgls_mwni(d,nh,nk,nx,map,w,niter,F)
eps=1e-3;
x=zeros(nx,1);
X=zeros(nk,1);
d=d(:);
y=[d;X];
dftop = 'no ';
if (nargin == 8) dftop = 'yes';end

if (dftop == 'yes')
    y0=sampling(X,nh,nk,nx,map,w,F);
else
    y0=sampling(X,nh,nk,nx,map,w);
end

y0(nh+1:nh+nk)=eps*X(1:nk);
rr=y-y0;

energy=rr'*rr;
e_old=energy;

if (dftop == 'yes')
    display('Using DFT');
    G=interpolate(rr,nh,nk,nx,map,w,F);
else
    display('Using FFT');
    G=interpolate(rr,nh,nk,nx,map,w);
end
    
G=G+eps*rr(nh+1:nh+nk);

S=G;
gammam=G'*G;
rms=0;

for iter=1:niter
    
    if (dftop == 'yes')
        ss=sampling(S,nh,nk,nx,map,w,F);
    else
        ss=sampling(S,nh,nk,nx,map,w);
    end
    
    ss(nh+1:nh+nk)=eps*S(1:nk);
    den=ss'*ss;
    %denp(iter)=den;plot(denp);figure(gcf);pause
    %if (abs(den)<1e-7) keyboard;end
    alpha=gammam/den;
    X=X+alpha*S;
    rr=rr-alpha*ss;
    if (dftop == 'yes')
        G=interpolate(rr,nh,nk,nx,map,w,F);
    else
        G=interpolate(rr,nh,nk,nx,map,w);
    end
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

X=X(:).*w(:);

if (dftop == 'yes')
    x=idft(F,X);
    x=x(1:nx);
    %update weights,
    XX=dft(F,x);

else
    x=ifft(X);
    x=x(1:nx);
    %update weights,
    XX=fft(x,nk);
end

w=abs(XX);w=smoothing(w,5);
%w=ones(size(w));

 %figure;plot(rmsv);

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
 
 display(sprintf('Repeated bins: removed %d ',nh-nh2));
 %ii=1:nh-1;find((map(ii+1)-map(ii))==0);
 return;
 
 function [F,w]=dftoperator(t,nz,wxi);
     % [F]=dftoperator(t,nz,nxi);
     % with t space axis (can be irregular)
     %      nz number of zeros to pad
     %      wxi optional weights proportional to the distance 
     %          between samples.
     %  Output; F is DFT operator.
     %          w is frequency axis.
     % Use this operator as follows:
     % DFT: X=F*x(:);
     % IDFT: x=F'*X(:);
     if (nargin < 3) wxi=ones(size(t));end
     if (nargin < 2) nz=0;end
     nt=length(t);
     dt=(t(end)-t(1))/(length(t)-1);
     
     % padding
     ntz=nt+nz;
     
     % generate a frequency axis that matches fft axis.
     f=(-ntz+2)/2:(ntz)/2;
     vv=fftshift(f);f=[vv(end) vv(1:end-1)];
     f=f/(dt*ntz);
     w=f*2*pi;
     %w=-lw/2:lw/2-1;
     %w=w/(lw-1)*2*pi;
     F=exp(-i*(w(:)*t(:).'))*diag(wxi);
     return;
 
 
 function [X]=dft(F,x)
    X=F*x(:);
    return

function [x]=idft(F,X)
    x=F'*X(:)/length(X);    
    return

function [FR]=reg_operator(t,w);
     % Uses a given frequency axis rather than calculating it from the
     % space axis. Use to reconstruct in a regular grid from a
     % spectrum calculated with a DFT operator.

     % [FR]=reg_operator(t,w);
     % Use this operator as follows:
     % DFT: X=FR*x(:);
     % IDFT: x=FR'*X(:);
     % In this example X has a frequency axis = w.
     nt=length(t);
     dt=(t(end)-t(1))/(length(t)-1);
          
     FR=exp(-i*(w(:)*t(:).'));
return;

function [data]=my_function(h)
    A1=1;
    A2=2;
    f1 = 0.05;
    f2 = 0.2;
    nh=length(h);
    data = A1*sin(2*pi*f1*h)+A2*sin(2*pi*f2*h);
        
return;

function plotspectrum(x,dx)
        plot(freqaxis(dx,length(x)),fftshift(abs(fft(x))));
return;
     
function [f]=freqaxis(dt,nt)
% function [f]=freqaxis(dt,nt)
f=(-nt+1)/2 - 1/2:(nt-1)/2 - 1/2;
f=f/(dt*nt);
return;

