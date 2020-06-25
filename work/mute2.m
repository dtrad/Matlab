function [vr]=mute2(vr,alfai,alfaf,ti,tf,dt)
global alfa t
% 
% input 
%  	vr : the original data, 
%     alfai : initial alfa
%     alfaf : final alfa 
%     ti : initial time
%     tf : final time 
% output 
% 		vr :  muted data 

% Daniel Trad-- 6-04-98

% Transform values of alfa_max alfa_min and t_max t_min to indexes
[nt nh]=size(vr),if nt < nh vr=vr.';end;
tti=round(ti/dt)+1;
ttf=round(tf/dt)+1;

indexmin=min(find(alfa>=alfai));
indexmax=max(find(alfa<=alfaf));

alfamute=indexmin:indexmax;
tt=tti:ttf;
vr(tt,alfamute)=zeros(size(vr(tt,alfamute)));
