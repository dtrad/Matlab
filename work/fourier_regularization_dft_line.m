function [log]=fourier_regularization_dft(doplot,method,firstfig)

% mwni interpolation for 2d
% irregular data testing 
ninterp=4;

% Parameters: 
%           doplot = 1: plot iteration by iteration
%           doplot = 2: normal plots on screen
%           doplot = 3: normal plots on screen and disk
%           method = 'reg', 'ireg'
%           firstfig = 1; use to preserve previous plots on screen.
% Daniel Trad; Veritas, March, 2006 (last May 4, 06)

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
%data=model(1:550,:);%clear model;
%h=xi;clear xi;
data=data(1:550,:);
t=t(1:550);
nt=length(t);dt=t(2)-t(1);


%figure, wigb(data,1,h,t);
%[data,h]=remove_traces(data,h);
%figure, wigb(data,1,h,t);

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
% regular axis (binning)
nh=length(h)
dx = dh/ninterp;
x=h(1):dx:h(end);
nx=length(x);
%%%%%%%%%%%%%%%%%%%%%%
% irregular axis for dft
index1=1:length(hi);
index2=1:1/ninterp:length(hi);
xi=interp1(index1,hi,index2,'spline');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate distance weights
wxi(1)=1;
for ix=2:nx-1
    % theoretical right weights. (doesn't work)
    %wxi(ix)=abs(xi(ix+1)-xi(ix-1));
    %wxi(ix)=wxi(ix)/(2*dx);
    % alternative weights (doesn't work).
    if (abs(xi(ix)-xi(ix-1)) < 1) 
        wxi(ix) = 1; % test different vales
    else
        wxi(ix) = 1;
    end    
end
wxi(nx)=1;
% Distance weights give bad results. Make them one for now.
wxi=ones(size(wxi));

% DFT operator
% 1) DFT for irregular axis.
[F,kx]=dftoperator(xi,nz,wxi);
% 2) DFT for regular axis.
[FR]=reg_operator(x,kx);

% Length for the kx axis 
nk=length(kx);


map=mapping(hi,xi);
model=upsample(data,map,nx);
display('size of model in t-x space');
size(model)

%zero padding and fft along time    
D=fft(data,2*nt);
fs=1/dt;nf=nt;
w=2*pi*(0:nf-1)*fs/nt;
DH=D(1:nf,:);     % data in f,x
MH=zeros(nf,nx);  % model in f,x
MHR=zeros(nf,nx); % model with regularization

weights=ones(1,nk); % CG weights.
% test the adjoint 
test=dotproduct_test([DH(1,:) zeros(1,nk)],fft(MH(1,:),nk),weights,nh,nk,nx,map,F);
display('test = ');test
niter=5;
[weights]=band_limit(weights,nk,ninterp);
log=zeros(nf,2);
for f=5:nf-50
    %For sinc interpolation leave weights=1
    % weights=ones(1,nk); % CG weights.

    [MH(f,:),weights,X,iter]=wtcgls_mwni(DH(f,:),nh,nk,nx,map,weights,niter,F);
    
    MHR(f,:)=idft(FR,X);
    % debugging plots to see effect of weigths
    if (doplot==1)
        subplot(311);plot(abs(fft(DH(f,:))));
        subplot(312);plot(abs(fft(MH(f,:))));
        subplot(313);plot(weights);figure(gcf);pause;
    end
    log(f,:)=[f iter];
end

log=log(5:nf-30,:);
M=duplic(MH);
MR=duplic(MHR);
model=ifft(M);model=model(1:nt,:);
modelr=ifft(MR);modelr=modelr(1:nt,:);
size(model)

ifig=firstfig;
if (doplot>1)
    figure(ifig+0),wigb(data,1,hi,t);
    if(doplot >2) print -dpng mwni_orig.png;end

    figure(ifig+1),wigb(model,ninterp,xi,t);
    if(doplot >2) print -dpng mwni_int.png;end
    FD=dftoperator(hi);
    subplot(111);figure(ifig+2),imagesc(fftshift(abs(dft2(FD,data))));
    if(doplot >2) print -dpng mwni_origspec.png;end
    clear FD;
    
    %subplot(111);figure(ifig+3),imagesc((abs(fft2(data))));
    
    FM=dftoperator(xi);
    figure(ifig+3),imagesc(fftshift(abs(dft2(FM,model))));
    if(doplot >2) print -dpng mwni_intspec.png;end
    clear FM;

    datap=downsample(model,map,nh);
    factor=(max(max(data-datap)))/(max(max(data)));
    figure(ifig+4),wigb(data-datap,factor,hi,t);
    if(doplot >2) print -dpng mwni_res.png;end

    figure(ifig+5),wigb(modelr,ninterp*0.5,x,t);
    if(doplot >2) print -dpng mwni_regint.png;end
end

return


function [d]=sampling(M,nh,nk,nx,map,weights,F);
    d=zeros(nh,1);
    weights=weights(:);
    M=M(:);
    M=M.*weights;
    m=idft(F,M);
    %m=ifft(M);
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
    M=dft(F,m)/nk;
    %M=fft(m,nk)/nk;
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
     return

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


function [x,w,X,iter]=wtcgls_mwni(d,nh,nk,nx,map,w,niter,F)
    eps=1e-5;
    x=zeros(nx,1);
    X=zeros(nk,1);
    d=d(:);
    y=[d;X];

    y0=sampling(X,nh,nk,nx,map,w,F);
    y0(nh+1:nh+nk)=eps*X(1:nk);
    rr=y-y0;

    energy=rr'*rr;
    e_old=energy;

    G=interpolate(rr,nh,nk,nx,map,w,F);
    G=G+eps*rr(nh+1:nh+nk);

    S=G;
    gammam=G'*G;

    for iter=1:niter
        ss=sampling(S,nh,nk,nx,map,w,F);
        ss(nh+1:nh+nk)=eps*S(1:nk);
        den=ss'*ss;
        %denp(iter)=den;plot(denp);figure(gcf);pause
        if (abs(den)<1e-14) 
            break;
        end
        alpha=gammam/den;
        X=X+alpha*S;
        rr=rr-alpha*ss;
        G=interpolate(rr,nh,nk,nx,map,w,F);
        G=G+eps*rr(nh+1:nh+nk);
        gamma=G'*G;
        if (abs(gammam)<1e-14) 
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
    x=idft(F,X);x=x(1:nx);

    %update weights,
    XX=dft(F,x);
    %XX=fft(x,nk);

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

function [test]=dotproduct_test(d,m,w,nh,nk,nx,map,F);
    D1=rand(size(d))+i*rand(size(d));
    M1=rand(size(m))+i*rand(size(m));
    D1=D1(:);
    M1=M1(:);
    weights=rand(size(w));
    [D2]=sampling(M1,nh,nk,nx,map,weights,F);
    D2(nh+1:nh+nk)=eps*M1(1:nk);

    [M2]=interpolate(D1,nh,nk,nx,map,weights,F);
    M2=M2+eps*D1(nh+1:nh+nk);

    dot1=D1(:)'*D2(:);
    dot2=M1(:)'*M2(:);

    test=abs(dot1/dot2);
    return

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

function [X]=dft2(F,x)
    [nk,nx2]=size(F);
    [nt,nx]=size(x);
    if (nt<nx)
        x=x.';
        [nt,nx]=size(x);
    end

    if (nx~=nx2)
        display('error: size of F and x are different');
                                                                return
    end

    xf=fft(x);nf=nt;
   
    X=zeros(nf,nk);
    for f=1:nf
        vtemp=xf(f,:);
        VT=F*vtemp(:);
        X(f,:)=VT(:).';
    end

    clear VT vtemp;

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

    function [model,x]=remove_traces(data,h)
        nh=length(h);
        dh = (h(end)-h(1))/(length(h)-1);
        
        model=zeros(size(data));
        
        ix=1;
        x(1)=h(1);
        model(:,1)=model(:,1);
        
        for ih=2:nh-1
            if (abs(h(ih+1)-h(ih-1)) > dh)
                
                x(ix)=h(ih);
                model(:,ix) =data(:,ih);
                ix=ix+1;
            end
        end
        x(ix)=h(nh);
        model(:,ix)=data(:,nh);
        nx=length(x);
        model=model(:,1:nx);
        
        return

        
        