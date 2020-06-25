function y=sinc(x)
ix=find(x==0);x(ix)=1;
y=sin(x)./x;
y(ix)=1;
return;
