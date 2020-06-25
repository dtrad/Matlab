function [cc]=reflectivity2(c,lx,nspikes)
if nargin<3 threshold=0.00001;end
dt=0.004;
cc=zeros(1,lx);
p=0.5;
for i=1:nspikes
   ind=random('unif',1,lx);
   nn=round(random('unif',1,2));
   ind=round(ind);sign=(-1)^nn;
   pp=random('unif',0,1);
   if pp>p 
      cc(ind)=sign*random('exp',0.01,1,1);
   else
      cc(ind)=sign*random('exp',0.08,1,1);
   end   
end;
indc=find(c~=0);
cc(indc)=c(indc);
Rx=autocorr(cc);
tt=0:length(cc)-1;tt=tt*dt;
tt2=0:length(Rx)-1;tt2=tt2*dt;tt2=tt2-max(tt2)/2;

figure
subplot(211);linesad(tt,cc);title('(a) Reflectivity');
subplot(212);linesad(tt2,Rx);title('(b) Autocorrelation');