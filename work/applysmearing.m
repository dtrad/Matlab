function [madj]=applysmearing(m,t,p,FF3)

dt=t(2)-t(1);
fs=1/dt;
nt=length(t);
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.

M=fft(m,nt);
[nf np]=size(M);
MADJ=zeros(nf/2,np);  

for f=2:nf/2-5,
  mm=M(f,:);
  adj=FF3(:,:,f)*mm(:);
  MADJ(f,:)=adj(:).';
end

%MADJ=zeropadm(MADJ,nt/2);
MADJ=duplic(MADJ);
madj=ifft(MADJ);

return;




