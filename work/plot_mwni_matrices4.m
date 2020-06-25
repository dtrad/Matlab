function plot_mwni_matrices3
%compare the system of equations for regular and irregular case    

tr=0:31;
bin = 0.5;
nt=length(tr);
baseWeight=99; % base weight for exponential weight (baseW to baseW + 1);
%create irregular axis
for i=1:nt
      ti(i) = tr(i) + 1*(rand() - 0.5);
end

% create regular axis for binning.
tb=0:bin:nt - 1;
nb = length(tb);

%amount of zero padding
nz = 0;

FR=dftoperator(tr,nz);
FI=dftoperator(ti,nz);
FB=dftoperator(tb,nz);


%generate some weight matrix as example 
wr=fftshift(baseWeight+exp(-(tr-nt/2).^2/100));
WR=diag(wr);

wi=fftshift(baseWeight+exp(-(ti-nt/2).^2/100));
WI=diag(wi);

wb=fftshift(baseWeight+exp(-(tb-nb/2).^2/100));
WB=diag(wb);


%Generate sampling matrices

% %%%%%%%%%%%%%%% gaps %%%%%%%%%%%%%%%%%%%%%
ii=floor(nt/4):floor(nt/2); % gap definition
T=eye(nt); % full sampled

T(ii,ii)=0; % gap;
TG=removeZeroRows(T); % sampling for data with gaps

% create data with the gap
tg=ti;
tg(ii)=0; 
tg=tg(find(tg ~= 0)); % remove zeros


%%%%%%%%%%%%%% decimation %%%%%%%%%%%%%%%%%%%%%
ii=2:2:nt; % decimation pattern
T=eye(nt);
T(ii,ii)=0;
TD=removeZeroRows(T); % sampling for data with decimation.

% create decimated data
td=ti;
td(ii)=0; 
td=td(find(td ~= 0));


%binning for gapped
tsg = floor(tg/bin)+1; % binned data.
tsg=tsg(find(tsg~=0)); % remove zeros
TBG=zeros(nb);
% build sampling with binning gapped axis 
for ii=1:length(tsg),TBG(tsg(ii),tsg(ii))=1;end
TBG=removeZeroRows(TBG);

%binning for decimated
tsd = floor(td/bin)+1; % binned data.
tsd = tsd( find (tsd ~= 0));

TBD=zeros(nb);
% build sampling with binning gapped axis 
for ii=1:length(tsd);TBD(tsd(ii),tsd(ii))=1;end;
TBD=removeZeroRows(TBD);


%figure(1);subplot(221);plot(tr,ti)
%figure(1);subplot(222);plot(tr,w);
%figure(1);subplot(223);imagesc(TD);
%figure(1);subplot(224);imagesc(TG);

LRD = TD*FR'*WR;
LRG = TG*FR'*WR;
LID = TD*FI'*WI;
LIG = TG*FI'*WI;

test = 0;
if ( test == 1) 
    % test, all bins are full.
    TBD = eye(nb);
    TBG = eye(nb);
    WB  = eye(size(WB));
end

LBD = TBD*FB'*WB;
LBG = TBG*FB'*WB;

ARD = abs(LRD'*LRD);
ARG = abs(LRG'*LRG);
AID = abs(LID'*LID);
AIG = abs(LIG'*LIG);
ABD = abs(LBD'*LBD);
ABG = abs(LBG'*LBG);

ARD=ARD/max(max(ARD));
ARG=ARG/max(max(ARG));
AID=AID/max(max(AID));
AIG=AIG/max(max(AIG));
ABD=ABD/max(max(ABD));
ABG=ABG/max(max(ABG));


CLIM=[0 1];
figure(2);
%colormap(jet);
subplot(321);imagesc(ARD,CLIM);colorbar;%xlabel('columm#');ylabel('row#')
subplot(322);imagesc(ARG,CLIM);colorbar;%xlabel('columm#');ylabel('row#')
subplot(323);imagesc(AID,CLIM);colorbar;%xlabel('columm#');ylabel('row#')
subplot(324);imagesc(AIG,CLIM);colorbar;%xlabel('columm#');ylabel('row#')
subplot(325);imagesc(ABD,CLIM);colorbar;%xlabel('columm#');ylabel('row#')
subplot(326);imagesc(ABG,CLIM);colorbar;%xlabel('columm#');ylabel('row#')

figure(3);
subplot(321);plot(sumDiag(ARD));xlabel('(a)');
subplot(322);plot(sumDiag(ARG));xlabel('(b)');
subplot(323);plot(sumDiag(AID));xlabel('(c)');
subplot(324);plot(sumDiag(AIG));xlabel('(d)');
subplot(325);plot(sumDiag(ABD));xlabel('(e)');
subplot(326);plot(sumDiag(ABG));xlabel('(f)');


figure(4);
% get axis from denser plot.
plot(eig(AIG));V=axis;


subplot(321);plot(eig(ARD));xlabel('(a)');axis(V);
subplot(322);plot(eig(ARG));xlabel('(b)');axis(V);
subplot(323);plot(eig(AID));xlabel('(c)');axis(V);
subplot(324);plot(eig(AIG));xlabel('(d)');axis(V);
subplot(325);plot(eig(ABD));xlabel('(e)');V(2)=V(2)/bin;axis(V);
subplot(326);plot(eig(ABG));xlabel('(f)');axis(V);
%plotAntiDiag(abs(LIG'*LIG),[0],'LIG',loc);
%plotAntiDiag(abs(LRG'*LRG),[0],'LRG'),loc;
%plotAntiDiag(abs(LBG'*LBG),[0],'LBG');


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
      
    function [a,la]=sumDiag(A)
    na=length(diag(A));
       
    for k=0:na-1,
        la(k+na) = length(diag(A,k));
        a(k+na)=sum(diag(A,k))/la(k+na);
        la(na-k) = length(diag(A,-k));
        a(na-k)=sum(diag(A,-k))/la(na-k);
    end;
    
            
    return
