function [vel,pr]=alfa2vel(alfa,p)
% [vel,pr]=alfa2vel(alfa,p)
% Map between alfa parameter and ray parameter (1/vel)
% Daniel Trad- UBC- 27-08-98
 
dp=p(2)-p(1)
dalfa=alfa(2)-alfa(1)
pr=t_2(dp,dalfa,alfa,-1);end
pr=sqrt(pr);
vel=1./(pr+1e-10);
