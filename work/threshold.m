function [T]=threshold(coef,perc)
% given a data set it finds the threshold to keep a certain percent
% or number of coefficients of coefficients 
if (perc<1)
  perc=perc*100;
  
  [N1 M1]=size(coef);
  total=N1*M1;
  kept=round((total)*(perc/100));
  if (kept>total) 
    kept=total;
  end
  message=sprintf('number of coefficients to keep=%d from %d or %f percent',...
		  kept,total,kept/total*100)
else
  kept=perc;
end

ss=sort((abs(coef(:))));T=ss(length(ss)-kept+1);

return;
