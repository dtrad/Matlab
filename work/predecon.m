function [xxp]=predecon(xx,xxm,delay,nn,h)
%function [xxp]=predecon(xx,delay,nn)
% Daniel Trad - Jun 4-98
xx=real(xx);
xxm=real(xxm);
axis1=[0 2.0 -20 1500];
[nt nh]=size(xxm);

if nargin<6 t=(0:nt-1)*0.004;end;
if nargin<5 h=(0:nh-1)*1; end; 

if max(size(delay))==1 delay(1:nh)=delay*ones(nh,1);end
if max(size(nn))==1 nn(1:nh)=nn*ones(nh,1);end

Rx=autocorr(xxm);%Rx(nt,1)=1.0;
%Rx(nt+11:nt+40,1)=zeros(size(Rx(nt+11:nt+40,1)));
nmax=nh;
%nmax=min(find(delay>=1))
%if (isempty(nmax)) nmax=nh;end

xxp=xx;

for hh=1:nmax;
   acoef=levinso1(Rx(nt:2*nt-1,hh),nn(hh),delay(hh));
   az=zeros(1,delay(hh)-1);
   ap=[1 az acoef(2:nn(hh))];
   ap=ap(:).';
   if hh==1 app=ap;end;
   temp=conv(xx(:,hh),ap);
	xxp(:,hh)=temp(1:nt);	
end;

xxtemp=xx./max(max(abs(xx)));
xxptemp=xxp./max(max(abs(xxp)));

figure
subplot(221);
wigb(real(xx),1,h,t);
title('original data')

%plotrace(xxtemp);
%axis(axis1);

subplot(222);
wigb(real(xxp),1,h,t);
title('predicted data')

%plotrace(xxptemp);
%axis(axis1);

subplot(223);
title('residuals');
wigb(real(xxtemp-xxptemp(1:nt,:)),0.5,h,t);

%plotrace(xxtemp-xxptemp(1:nt,:));
%axis(axis1);

subplot(224);
plot(Rx(nt:nt+nn(1)+delay,1));
title('autocorrelation function for trace #1 in tau-p domain')

%subplot(224);plotrace(Rx(nt:nt+nn+delay,:));title('autocorrelation function')

figure;
subplot(211);
wigb(detrend(real(Rx(nt:2*nt-1,:))),1,h,t);
%plotrace(detrend(Rx(nt:2*nt-1,:)));
title('autocorrelation function before pd')
VV=axis;
Rx=autocorr(xxp);nt=max(size(xxp));Rx(nt,1)
subplot(212);
wigb(detrend(real(Rx(nt:2*nt-1,:))),1,h,t);
%plotrace(detrend(Rx(nt:2*nt-1,:)));
title('autocorrelation function after pd')
axis(VV);
figure,
plot(app);title('filter coeff');

