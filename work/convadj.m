% test for conv/correlation adjoint pair
N=256;
t=1:256;
r=zeros(1,N);
r(1,100)=1;
w=ricker;
lf=length(w)
lfh=(lf-1)/2
x=conv(r,w);
x2=x(lfh+1:N+lfh);
x=x(1:N);

subplot(211),plot(t,r,t,x,t,x2);


 wn=w(end:-1:1);
 y=conv(x,wn);
 y2=y(lf:end);
 %plot(t,r,t,x,t,y2)
 subplot(212);plot(t,r,t,x,t,y2/10)
 
 
 figure(gcf);