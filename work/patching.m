function dout=patching(d,dout,f1,n1,f2,n2)
% function dout=window(d,f1,n1,f2,n2)

dout(f1:f1+n1-1,f2:f2+n2-1)=d;

return