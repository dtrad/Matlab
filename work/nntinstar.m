% 
% Example: Hebbian learning

V= [ 0.4082; 0.8165; -0.4082]


P=randn(3,200)*0.1+V*ones(1,200);

P=normc(P);

W=[-0.8133 0.1474 -0.5628]


A0=satlin(W*V)

lr=0.05;

for q=1:200
  p = P(:,q);
  a = satlin(W*p);
  dW= learnis(W,p,a,lr);
  W=W+dW;
end

W
V.'

A0
A=satlin(W*V)

