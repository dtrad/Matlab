function [p,dq]=radonaxis2(np,vmin,vmax,h,curve)
if (nargin < 5) curve='PRT';end

fmax=70;

dh=h(2)-h(1);
if (curve=='PRT')
   %pmin=1/(vmax^2);
   pmin=-2e-8
   dq=abs(0.8/(fmax*(max(h)-min(h))^2))
   pmax=pmin+(np-1)*dq
   p=pmin:dq:pmax;
   %pmaxt=1/(2*max(h)*fmax*dh)+pmin
elseif (curve=='HRT')
   pmin=1/(vmax^2);
   pmax=1/(vmin^2);
   dq=(pmax-pmin)/(np-1);
   p=pmin:dq:pmax;
elseif (curve=='LRT')
   pmin=0;
   pmax=1/(vmin);
   dq=(pmax-pmin)/(np-1);
   p=pmin:dq:pmax;   
end



