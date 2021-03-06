function [u]=backward2(V,UH,w,h1,nt,Power,WV,WU0,WU1,W,alfa)
% Backward Transform: u~=L*v or u~=F*WV.v
% u~ are the recovered data
% 		V w-p data
%   	UH w-x data
% output 
%   	u t-x data
%		JD data misfit

% Daniel Trad-- 6-04-98

h=h1;
WU=WU0;
JD=0;

for f=1:nt/2;
   F=exp(i*w(f)*(alfa*(h.^Power)));
   FH=F';
   ut=(FH*WV)*(V(f,:)).';
   if(f==1) UT=ut;
   else UT=[UT,ut];
   end,
   %utemp=shrinktr(ut,HH);utemp=utemp(:); 
   %if (f<=60) JD=JD+(W*(utemp-UH(f,:).'))'*WU*(W*(utemp-UH(f,:).'));end;
end;
UT=UT.';
UTD=duplic(UT);
u=ifft(UTD);

