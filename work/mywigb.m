function mywigb(R)
[NF,NH]=size(R);
figure,wigb(ifft(duplic(R(1:NF/2,:))));

