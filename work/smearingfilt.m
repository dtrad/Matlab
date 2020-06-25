function [F3,FF3]=smearingfilt(t,h,p,rtmethod,maxindex)
% [F3,FF3]=smearingfilt(t,h,p,rtmethod)
% input the three axes for the RT operator

if nargin<4|isempty(rtmethod) rtmethod='PRT';end  


nt=length(t);
np=length(p);
nh=length(h);



dt=t(2)-t(1);
fs=1/dt;  % sample freq.



w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.
nf=length(w);
if nargin<5|isempty(maxindex) maxindex=nf;end
if maxindex>nf;maxindex=nf;end

FF3=zeros(np,np,maxindex);
F3=zeros(nh,np,maxindex);



for f=2:maxindex,
   [FH,F]=radonmat2(w(f),h,p,rtmethod);
   if (nargout>1) 
     FF3(:,:,f)=FH*F;
     display('calculating  smearing filter as well');
   end  
   F3(:,:,f)=F;
   
     
end


return;