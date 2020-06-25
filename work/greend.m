function [G0]=greend(xs,zs,kxg,zg,q,rho)
% Green's function for the wave equation in a halfspace
% [G0]=greend(xs,zs,kxg,zg,q,rho)


lw=max(size(q));
kxg=kxg(:).';

cte=-rho/2/sqrt(2*pi)/i;
G0=cte*exp(-i*(ones(lw,1)*kxg*xs+q.*(ones(lw,1)*abs(zg-zs))));
G0=G0./q;