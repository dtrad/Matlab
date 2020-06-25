function [madj]=applyadj(d,t,p,F3)

dt=t(2)-t(1);
fs=1/dt;
nt=length(t);
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.

D=fft(d,nt);
[nf nh]=size(D);
np=length(p);
M=zeros(nf/2,np);  

for f=2:nf/2-5,
  dd=D(f,:);
  mvector=F3(:,:,f)'*dd(:);
  M(f,:)=mvector(:).';
end

%MADJ=zeropadm(MADJ,nt/2);
M=duplic1(M);
madj=ifft(M);

return;




