%function [wav]=wavelet2d(sx,st,tau,nt,nx,dt,dx)

sx=600;
st=0.005;
tau=1;
nt=512;
nx=40;
dt=0.004;
dx=25;

x=-nx/2:nx/2;
x=x*dx;x=x(:).';
t=-nt/2:nt/2;
t=t*dt;t=t(:);

gridt=t*ones(size(x));
gridx=ones(size(t))*x;

%tmp=pow(((t-sqrt((x/sx)*(x/sx)+tau*tau))/st),2);
%tmp=((gridt-sqrt((gridx/sx).*(gridx/sx)+tau*tau))/st);
%wav=((1-tmp).*exp(-tmp/2))/(sx*st);


for it=1:nt
   for ix=1:nx
     tmp(it,ix)=((t(it)-sqrt((x(ix)/sx).*(x(ix)/sx)+tau*tau))/st);
     wav(it,ix)=((1-tmp(it,ix)).*exp(-tmp(it,ix)/2))/(sx*st);
   end
end


wigb(wav);
