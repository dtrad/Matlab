function [u]=backward(V,UH,w,h1,nt,Power,WV,WU0,WU1,W,alfa,HH)
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
end;
UT=UT.';
UTD=duplic(UT);
u=ifft(UTD);

