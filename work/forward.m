function [V]=Forward(UH)
% Forward Transform: v=Lu or v=FWU.u
% u are the original data, 
% input u t-x data
% output 
% V  w-p data 

% Daniel Trad-- 6-04-98
global w dte h0 nt Power WV WU0 alfa

h=h0;
WU=WU0;

for f=1:nt/2;
   F=exp(i*w(f)*(alfa*(h.^Power)));
   v=(F*WU)*(UH(f,:)).';
   if(f==1) V=v;
   else V=[V,v];
   end,   
end;
