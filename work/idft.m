function x=idft(X,t,w);
lx=length(t);
lw=length(w);
X=X(:);

dt(1)=0;
dt(2:lx)=t(2:lx)-t(1:lx-1);

F=exp(-i*(t(:)*w(:).'));

x=F*X;