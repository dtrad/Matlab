function simage(D,x,t)
% SIMAGE    SIMAGE(D,x,t) plots a seismic image.
%           The image is given as a matrix D and the axis
%           are plotted according to the vectors t and x.
%
%
%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%
%  sacchi@phys.ualberta.ca
%

nt=max(size(D));
nx=min(size(D));


if (nargin < 1 | nargin > 3)
 error('Wrong number of input parameters.')
  end     

if (nargin < 3)
 x=1:1:nx;
  end     

if (nargin < 2)   % assume 4msec
 dt=0.004;
  t=0:dt:(nt-1)*dt;
   end     

Dmin=min(min(D)); Dmax=max(max(D));
Dn=(D-Dmin)*64/(Dmax-Dmin);

image(x,t,Dn);
