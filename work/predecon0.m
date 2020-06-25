function [xxp]=predecon(xx,xxm,delay,nn)
%function [xxp]=predecon(xx,delay,nn)
% Daniel Trad - Jun 4-98
xx=real(xx);
xxm=real(xxm);
axis1=[0 2.0 -20 1500];
[nt nh]=size(xxm);
if max(size(delay))==1 delay(1:nh)=delay*ones(nh,1);end
Rx=autocorr(xxm);%Rx(nt,1)=1.0;
%Rx(nt+11:nt+40,1)=zeros(size(Rx(nt+11:nt+40,1)));

nmax=min(find(delay==0));

xxp=xx;

for hh=1:nmax;
   acoef=levinso1(Rx(nt:2*nt-1,hh),nn,delay(hh));
   az=zeros(1,delay(hh)-1);
   az2=zeros(1,delay(1)-delay(hh));
   if delay(hh)>0
      ap=[1 az acoef(2:nn) az2];
      ap=ap(:).';
      if hh==1 app=ap;end;
   elseif delay(hh)==0;
      az1=acoef(2:nn).*0;
      ap=[1 az az1 az2];
   end;
   temp=conv(xx(:,hh),ap);
	xxp(:,hh)=temp(1:nt);	
end;

xxtemp=xx./max(max(abs(xx)));
xxptemp=xxp./max(max(abs(xxp)));

figure
subplot(221);plotrace(xxtemp);title('original data')
axis(axis1);
subplot(222);plotrace(xxptemp);title('predicted data')
axis(axis1);
subplot(223);plotrace(xxtemp-xxptemp(1:nt,:));title('residuals')
axis(axis1);
subplot(224);plot(Rx(nt:nt+nn+delay,1));title('autocorrelation function')
%subplot(224);plotrace(Rx(nt:nt+nn+delay,:));title('autocorrelation function')
figure;
subplot(211);plotrace(detrend(Rx(nt:2*nt-1,:)));title('autocorrelation function before pd')
VV=axis;
Rx=autocorr(xxp);nt=max(size(xxp));Rx(nt,1)
subplot(212);plotrace(detrend(Rx(nt:2*nt-1,:)));title('autocorrelation function after pd')
axis(VV);
figure,
plot(app);title('filter coeff');

