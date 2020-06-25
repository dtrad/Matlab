function [unmo]=rnmo3(u,vel,h0,dt,dir)
% NMO correction
% input 
%  	u  : the original data, 
%     vel: velocity 
%		h0 : offset
%     dt : sample time
% output 
% 		unmo :  output data
%
% Daniel Trad-- 6-06-98

if nargin < 4 dt=0.004;end;
if nargin<5 dir=+1;end

[nt nh]=size(u);
unmo=zeros(size(u));

for hh=1:nh;
   h=h0(hh);
	for tt0=1:nt;
   	t0=(tt0-1)*dt;
   	t=sqrt(t0.^2+(h/vel).^2);
      tti=t./dt+1;
		   
      AMP2=tti-floor(tti);
      %AMP2=0;
      AMP1=1-AMP2;
      
      tti=floor(tti);
      %tti=round(tti);

      dti=tti-tt0;
      
      if dir >0 & (tt0+dir*dti+1) <= nt
         unmo(tt0,hh)=AMP1*u(tt0+dir*dti,hh)+AMP2*u(tt0+dir*dti+1,hh);
      end
      
      if  dir<0 & (tt0+dir*dti) > 0 & (tt0+dir*dti+2) <= nt
         unmo(tt0,hh)=AMP1*u(tt0+dir*(dti-1),hh)+AMP2*u(tt0+dir*(dti-1)+1,hh);
      end

      
	end
end
