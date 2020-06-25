x=ones(30,1);
for i=1:5,
  x(i)=1-exp(-(i-1)+.3);
  x(30-i+1)=x(i);
end,
