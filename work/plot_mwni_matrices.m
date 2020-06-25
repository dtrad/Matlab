function plot_mwni_matrices
    
t=0:31;
nt=length(t);
irreg = 1;

if (irreg == 1)
    for i=1:nt
        t(i) = t(i) + 1*(rand() - 0.5);
    end
end
    
nz = 0;
F=dftoperator(t,nz);
w=fftshift(1+exp(-(t-nt/2).^2/100));
W=diag(w);
I=eye(nt);

%gaps
T=I;for ii=floor(nt/4):floor(nt/2);T(ii,ii)=0;end
%decimation
%T=I;for ii=1:2:nt;T(ii,ii)=0;end

TT=[];
for i=1:length(T(:,1))
    if (T(i,i)~= 0) 
        TT=[TT;T(i,:)];
    end;
end

T=TT;

figure(1);subplot(311);plot(t)
figure(1);subplot(312);plot(w);
figure(1);subplot(313);plot(diag(T));

FF=F*F';
FTTF=F*T'*T*F';
WFFW=W'*F*F'*W;
WFTTFW=W'*F*T'*T*F'*W;

figure(2);
colormap(jet);
subplot(221);imagesc(abs(FF));colorbar;xlabel('(a)');
subplot(222);imagesc(abs(FTTF));colorbar;xlabel('(b)');
subplot(223);imagesc(abs(WFFW));colorbar;xlabel('(c)');
subplot(224);imagesc(abs(WFTTFW));colorbar;xlabel('(d)');

figure(3);
subplot(221);plot(diag((abs(FF))));xlabel('(a)');
subplot(222);plot(diag(abs(FTTF)));xlabel('(b)');
subplot(223);plot(diag(abs(WFFW)));xlabel('(c)');
subplot(224);plot(diag(abs(WFTTFW)));xlabel('(d)');

figure(3);
subplot(221);plot(antidiag((abs(FF))));xlabel('(a)');
subplot(222);plot(antidiag(abs(FTTF)));xlabel('(b)');
subplot(223);plot(antidiag(abs(WFFW)));xlabel('(c)');
subplot(224);plot(antidiag(abs(WFTTFW)));xlabel('(d)');
keyboard

return;
figure(4);
colormap(summer);
subplot(221);mesh(abs(F*F'));colorbar;xlabel('(a)');
subplot(222);mesh(abs(F*T*T*F'));colorbar;colorbar;xlabel('(b)');
subplot(223);mesh(abs(W*F*F'*W));colorbar;colorbar;xlabel('(c)');
subplot(224);mesh(abs(W*F*T*T*F'*W));colorbar;colorbar;xlabel('(d)');

return;
function [F,w]=dftoperator(t,nz,wxi)
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
    function [B]=antidiag(A)
        B=diag(A(end:-1:1,:));
    return;
        