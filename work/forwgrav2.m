function [F]=forwgrav2(z,N,M,delx,delz)
n0=(N-1)/2+1;
m0=(M)/2;
for i=2:M+1,
  for j=1:N,
    x=(-n0+j+m0-(i-1))*delx;
    x=x^2;
    top=x+z^2+z*delz+delz^2/4.0;
    bot=x+z^2-z*delz+delz^2/4.0;
    F(i-1,j)=4.68*log(top/bot);
  end,
end,
