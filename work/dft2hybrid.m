function [X]=dft2hybrid(x,h)
    F=dftoperator(h);
    size(F)
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
    xfh=xf(1:nt,:);
    X=zeros(nf,nk);
    for f=1:nf
        vtemp=xfh(f,:);
        VT=F*vtemp(:);
        X(f,:)=VT(:).';
    end

    clear VT vtemp;

    return
    
function [F,w]=dftoperator(t);
% [F]=dft(x,t);
% Use this operator as follows:
% DFT: X=F*x(:);
% IDFT: x=F'*X(:);
nt=length(t);
dt=(t(end)-t(1))/(length(t)-1);
f=(-nt+2)/2:(nt)/2;
vv=fftshift(f);f=[vv(end) vv(1:end-1)]
%f=[0:nt/2 (nt-1)/2:-1:1]; 
f=f/(dt*nt);
w=f*2*pi;
%w=-lw/2:lw/2-1;
%w=w/(lw-1)*2*pi;
F=exp(-i*(w(:)*t(:).'));
return;