% 
% Example: Hebbian learning


P=[1 0 0 0 1;
   0 0 0 0 0;
   0 0 1 1 0;
   0 0 0 0 1;
   1 0 0 0 1];

index=floor(rand(1,50)*4)+1;

P=P(:,index);

W=eye(5);
b=-0.5*ones(5,1);

figure(1);
hintonwb(W,b);

lr=1;
dr=0.01;
for q=1:50
  p = P(:,q);
  a = hardlim(W*p,b);
  dW= learnhd(W,p,a,lr,dr);
  W=W+dW;
end

figure(2)
hintonwb(W,b);