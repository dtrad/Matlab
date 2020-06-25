function [taper]=mytaper(ntaper)

for k = 1:ntaper 
  s = sin(k*pi/(2*ntaper));
  %taper(k) = s*s;
  %s = exp(-exp(k*pi/(2*ntaper)));
  taper(k)= s^100;%*s*s*s*s*s;
end


