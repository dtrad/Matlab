function w = ricker(f,dt)
% RICKER(f,dt)
% Returns a Ricker wavelet of
% central freq. f.   
% dt is the sampling interval


% Mauricio D. Sacchi
% Jan/1997.
% sacchi@geop.ubc.ca

% f << 1/(2dt)....
if (nargin<2) dt=0.004;end
if (nargin<1) f=25;end

nw=4./f/dt;
nw=2*floor(nw/2)+1

nc=floor(nw/2)
i=1:nw;
alpha=(nc-i+1).*f*dt;
beta=alpha.^2;
w=(1.-beta.*2).*exp(-beta);

w=w(:);