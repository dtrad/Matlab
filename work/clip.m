function  Dc  = clip(D,perc_low, perc_upper);
% CLIP     CLIP(D,perc_low,perc_upper) clips the data.
%          perc_low is the percentage of clipping for
%          negative amplitues and perc_upper is the percentage
%          of clipping for the positive amplitudes.
%
%          Note: the clipping should be done after looking
%          the histogram to decide where to clip.
%
%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%
%  sacchi@phys.ualberta.ca
%


	if (nargin < 1 | nargin > 3)
         error('Wrong number of input parameters in CLIP.')
	  end

	Dmax = max(max(D));
	Dmin = min(min(D));

	clip_low = Dmin-Dmin*perc_low/100.;
	clip_upper  = Dmax-Dmax*perc_upper/100.;

	[i]=find(D>=clip_upper);
	D(i) = clip_upper;

	[k]=find(D<=clip_low);
	D(k) = clip_low;

        Dc = D;
