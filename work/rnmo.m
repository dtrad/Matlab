function [unmo]=rnmo(u,nto)
% 
% input 
%     ti : initial time for trace 1,
%     tf : final time for last trace
%  	ur : the original data, 
% output 
% 		ur :  muted data 

% Daniel Trad-- 6-04-98
global w dte p h t np nh nt Power WV WU W dh h_near ff tor dp dt
unmo=zeros(size(u));

for hh=1:nh;
   tti(hh)=max(find(u(:,hh)>0.01));
   dti(hh)=tti(hh)-tti(1);
   for tt=1:nto-dti(hh)
   unmo(tt,hh)=u(tt+dti(hh),hh);
	end
end
tti