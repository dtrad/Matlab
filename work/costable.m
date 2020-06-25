function [yt yc]=costable(xx,a,b)
%vector to compute sin/cos
%a sin
%b cos

n=2881;
h=1440;
h1=1441;
s=zeros(n,1);
c=zeros(n,1);
d=pi/(8*180);

for i=1:h
  %[h-i h+i]
  x=(i-1)*d;
  s(h+i)=-sin(x);
  s(h-i+1)=sin(x);
  c(h+i)=-cos(x);
  c(h-i+1)=c(h+i);
end

if (a)
  yt=sin(xx)';
  yc=mytable(xx,d,s,n,h);
end

if (b)
  yt=cos(xx)';
  yc=mytable(xx,d,c,n,h);
end




function [y]=mytable(x,d,s,n,h)
y=0;
m=mod(x,2*pi);
i=round(m./d)+1;
i(find(i<1))
i(find(i>n))
y=s(i);
return;

