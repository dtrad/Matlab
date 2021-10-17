function [y]=nmo(x,dt,h,v)
% function [y]=nmo(x,dt,h,v)
% Daniel Trad - 
y=zeros(size(x));
nt=length(x);
for it0=1:nt
  t0=it0*dt;
  t=sqrt(t0^2+h^2/v^2)
  it=round(t/dt+1.5);
  if (it<nt) y(it0)=x(it);end
  % if synthetic y[it]=x(it0)
end
return;

function vel=interpolation(tnmo, vnmo)
for it=0:nt
  a=tnmo(t)/dt-it
  vel(it)=(1-a)*vnmo(t)+a*vnmo(t+1)
end
%size of vel = ntx1

vel=interpolation(tnmo,vnmo)

h=[hmin:dh:hmax];
for ih=1:nh
  x=data(ih,:)
  y=nmo(x,dt,h[ih],vel)
  dataout(ih,:)=y(:)
end
