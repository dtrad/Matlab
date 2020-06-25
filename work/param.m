function [alfa,dalfa]=param(optionp,p0,pmax,np)
% Computes the parameter alfa for the three possible cases:
% Input:
%   optionp:'Hampson ', 'Yilmaz  ', 'linear  '
%   p0   (usually 0)
%   pmax (1/vel max)
%   np    
% ouput
%    alfa  parameter for velocities
%    dalfa alfa interval
% Daniel Trad UBC, 6-07-98
dp = (pmax-p0)/(np-1);
p = p0:dp:pmax; 

if optionp=='Yilmaz '  
   alfamin=p0.^2;alfamax=pmax.^2;
   dalfa=(alfamax-alfamin)/(np-1);
   alfa=alfamin:dalfa:alfamax;alfa=alfa(:);
end
if optionp=='Hampson' alfa=p;alfa=alfa(:);dalfa=dp;end
if optionp=='linear ' alfa=p;alfa=alfa(:);dalfa=dp;end
