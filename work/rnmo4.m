function [unmo]=rnmo(u,vel,t0,dt,dh,h_near)
% 
% input 
%     t0 : time at zero offset
%     vel: velocity 
%  	u  : the original data, 
% output 
% 		unmo :  output data

% Daniel Trad-- 6-06-98
if nargin < 6 h_near=0;end,
if nargin < 5 dh=50;end,
if nargin < 4 dt=0.004;end,

[nt nh]=size(u);if nt<nh u=u.';end;
[nt nh]=size(u);

unmo=zeros(size(u));

tt0=round(t0/dt)+1;
for hh=1:nh;
   h=dh*(hh-1)+h_near;
   t=sqrt(t0.^2+(h/vel).^2);
   tti(hh)=round(t./dt)+1;
   dti(hh)=tti(hh)-tt0;
   for tt=1:nt-dti(hh)
   unmo(tt,hh)=u(tt+dti(hh),hh);
	end
end

