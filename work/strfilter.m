function [df]=strfilter(x,dt,h,v)
% function [y]=nmo(x,dt,h,v)
% Daniel Trad - 
df=zeros(size(x));
nt=length(x);
%keyboard
for it0=1:nt
  t0=it0*dt;
  t=sqrt(t0^2+h^2/v^2);
  dtnmo=t-t0;
  df(it0)=dtnmo/t0;
end

return;
