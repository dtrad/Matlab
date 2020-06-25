function [x]=integrate2fft(x)
n=length(x);
m=2;
while ( m < n) m=m*2;
end
y=padzeros(x,m);
dt=0.004;
w=freqaxisCalc(dt,m);

Y=fft(y);
for i=2:m
  Y(i)=Y(i)/(w(i)*w(i));
end
x=real(ifft(Y));
x=x(1:n)*m*4;
x=x-mean(x);

return;

  
