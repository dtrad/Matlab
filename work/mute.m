function [ur]=mute(ur,ti,tf)
% 
% input 
%     ti : initial time for trace 1,
%     tf : final time for last trace
%  	ur : the original data, 
% output 
% 		ur :  muted data 

% Daniel Trad-- 6-04-98
global w dte p h t np nh nt Power WV WU W dh h_near ff tor dp dt
ttf0=round(tf/dt)+1;
tti0=round(ti/dt)+1;

for hh=1:nh;
   tti(hh)=round(hh*(ttf0-tti0)/nh+tti0);
   ur(1:tti(hh),hh)=zeros(size(ur(1:tti(hh),hh)));
end
