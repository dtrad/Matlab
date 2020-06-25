% 
% Example: Hebbian learning
clear


P1=[1 0];
P2=[5 4;
    -2 1;
    3 -6];

W=zeros(3,1);


lr=1;

for q=1:2
  p1 = P1(:,q);
  p2 = P2(:,q);
  a = purelin(W*p1+p2)
  dW= learnh(p1,a,lr);
  W=W+dW
end

W



A=purelin(W*1+0);

