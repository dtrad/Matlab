function [fp]=firstpeak(xx,delay)
% function [xxp]=firstpeak(xx,delay)
% Daniel Trad - Jun 4-98
xx=real(xx);
axis1=[0 2.0 -20 2000];
[nt nh]=size(xx);
xxfp=zeros(size(xx));
for hh=1:nh;
   fp(hh)=find(max(abs(xx(delay:nt,hh))));
   xxfp(fp(hh))=xx(fp(hh),hh);
end;

xxtemp=xx./max(max(abs(xx)));
xxfptemp=xxfp./max(max(abs(xx)));

figure
subplot(211);plotrace(xxtemp);title('original data')
axis(axis1);
subplot(212);plotrace(xxfptemp);title('first peak')
axis(axis1);

