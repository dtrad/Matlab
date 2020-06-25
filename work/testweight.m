function [x] = testweight(threshold)
%clear;
%threshold=5;
nx=10;
x=[0:nx-1];

tolx=3; 
 
for i=1:nx
  expo1(i) = threshold * (abs(x(i)) - threshold);
  if (expo1(i) > tolx) 
    w(i)=1;
  else
    w(i)= 1 - exp(-exp(expo1(i)));
  end
  %[i x(i) expo1 w(i)] 

end  
%return
ii=1:nx;
constant=ones(1,nx);
plot(ii,x,ii,w,ii,constant*threshold,ii,expo1,ii,exp(expo1))
