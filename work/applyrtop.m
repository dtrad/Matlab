function [d]=applyrtop(m,t,h,F3)

[N1 N2 N3]=size(F3);

dt=t(2)-t(1);
fs=1/dt;
nt=length(t);
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.

M=fft(m,nt);
[nf np]=size(M);
nh=length(h);
D=zeros(nf/2,nh);  

for f=2:N3,
  mm=M(f,:);
  dvector=F3(:,:,f)*mm(:);
  D(f,:)=dvector(:).';
end

%MADJ=zeropadm(MADJ,nt/2);
D2=duplic1(D);
d=ifft(D2);

return;




