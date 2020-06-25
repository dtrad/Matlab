function plot_mwni_matrices2
%compare the system of equations for regular and irregular case    
tr=0:31;
nt=length(tr);
%create irregular axis
for i=1:nt
      ti(i) = tr(i) + 1*(rand() - 0.5);
end
%amount of zero padding    
nz = 0;
FR=dftoperator(tr,nz);
FI=dftoperator(ti,nz);

%generate some weight matrix as example (for both cases)
w=fftshift(1+exp(-(tr-nt/2).^2/100));
W=diag(w);

%Generate sampling matrices
I=eye(nt);
%gaps
T=I;for ii=floor(nt/4):floor(nt/2);T(ii,ii)=0;end
TD=removeZeroRows(T);

%decimation
T=I;for ii=1:2:nt;T(ii,ii)=0;end
TG=removeZeroRows(T);

figure(1);subplot(221);plot(tr,ti)
figure(1);subplot(222);plot(tr,w);
figure(1);subplot(223);imagesc(TD);
figure(1);subplot(224);imagesc(TG);

LRD = TD*FR'*W;
LRG = TG*FR'*W;
LID = TD*FI'*W;
LIG = TG*FI'*W;


figure(2);
%colormap(jet);
subplot(221);imagesc(abs(LRD'*LRD));colorbar;xlabel('(a)');
subplot(222);imagesc(abs(LID'*LID));colorbar;xlabel('(b)');
subplot(223);imagesc(abs(LRG'*LRG));colorbar;xlabel('(c)');
subplot(224);imagesc(abs(LIG'*LIG));colorbar;xlabel('(d)');

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
    
    function [B]=removeZeroRows(A)
    nz=find(sum(A,2));
    B=A(nz,:);
    return;
      
