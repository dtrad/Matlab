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

function 
