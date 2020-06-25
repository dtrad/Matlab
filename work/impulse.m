function [xx]=impulse(S,it,ip,NT,NX)
[M N]=size(S);
ipt=(ip-1)*NT+it;
m=sparse(ipt,1,1,N,1);
x=S*m;
xx=reshape(x,NT,NX);
wigb(full(xx));

