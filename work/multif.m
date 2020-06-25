function [t]=multif(t0,h,v,xm)

beta=33.3;
Rn=1840;
Rnip=184;

t=sqrt( (t0 + 2*xm*sin(beta)/v).^2 + 2*t0*cos(beta).^2/v*(xm.^2/Rn + h.^2/Rnip));

