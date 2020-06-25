function [p]=radonaxis2(np,vmin,vmax,h)
fmax=125;
dh=h(2)-h(1)
pmin=1/(vmax^2);
dq=0.8/(fmax*(max(h)^2-min(h)^2))
pmax=pmin+(np-1)*dq
p=pmin:dq:pmax;
pmaxt=1/(2*max(h)*fmax*dh)+pmin
