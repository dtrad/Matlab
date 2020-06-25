function [vr]=rmsvel(v,t)
% Rms velocities from interval vel v  and time t
% [vr]=rmsvel(v,t)
% where v is av ector with interval velocities and vr is vector
% with rms vel.
% t contains the total two way travel time from surface to the 
% interface (it is not the thickness but depth in time)
% Daniel Trad - UBC

nt=length(t);

nv=length(v);

% if (nt~=nv) display('nt must be equal to nv');

dt(1)=t(1);

for i=2:nt
 dt(i)=t(i)-t(i-1)
end

dt

sum=v(1)^2*dt(1);
ttot=dt(1);
vr(1)=v(1);
for i=2:nv
 ttot=ttot+dt(i);
 sum=sum+v(i)^2*dt(i);
 vr(i)=sqrt(sum/ttot);
end
