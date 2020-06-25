function [m]=applydesmearing(madj,t,p,FF3,eps1)

dt=t(2)-t(1);
fs=1/dt;
nt=length(t);
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.

M=fft(madj,nt);
[nf np]=size(M);
MADJ=zeros(nf/2,np);  

for f=2:nf/2-5,
  mm=M(f,:);
  adj=inv(FF3(:,:,f)+eps1*eye(np))*mm(:);
  MADJ(f,:)=adj(:).';
end

%MADJ=zeropadm(MADJ,nt/2);
MADJ=duplic1(MADJ);
m=ifft(MADJ);

return;




