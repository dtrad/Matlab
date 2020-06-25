% Program example to study 2D fft
NR=8;NC=8;
a=zeros(NR,NC);
for ii=1:NR,a(ii,ii)=ii;end;
af=fft2(a);
Q=1;
while (Q~=0)
Q=input('quadrant?');
if Q==1
   af(NR/2+1:NR,1:NC/2)
elseif Q==2
   af(1:NR/2,1:NC/2)
elseif Q==3
   af(1:NR/2,NC/2+1:NC)
elseif Q==4
   af(NR/2+1:NR,NC/2+1:NC)
end
end
