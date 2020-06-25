function y = ThresBand(x,dlevs,rlevs,par,N,nvar);
% This function removes the bands specified in rlevs
% INPUT
%            x     -    cell array with the curvelet decomposition
%            dlevs -    levels of the curvelet decomposition
%            rlevs -    cell array with levels to be zeroed
% Output
%            y     -     cell array with zeroed bands
  
   if length(dlevs)~=length(rlevs)-1,
 
     error('The size of vector rlevs should equal size of dlevs');
     return;
    
   end

  y  = x;

  for iband = 1 : length(dlevs)
    for iiband = 1 : 2^dlevs(iband)
% $$$       if (size(rlevs{iband})+1~= 2^(dlevs(iband)-1)),
% $$$  	error('not the right number of angles in rlevs')
% $$$  	return;
%      end
      if rlevs{iband+1}{iiband}==0,
	y{iband+1}{iiband} = zeros(size(y{iband+1}{iiband}));
      end
% $$$       if  rlevs{iband+1}{iiband}<0,
% $$$ 	y{iband+1}{iiband} = threshold_band(x{iband+1}{iiband},dlevs,rlevs{iband+1}{iiband});
%      end
    end
  end
  
  if nargin >3
    x = y;
    
    fac         =    1;
    
    pdfb_thr    =    fac * par * 3 *  nvar;

    [vcx,s]     =    pdfb2vec(x);
    
    fssize      =    prod(s(end, :));	% finest scale size
    
    pdfb_thr(end-fssize+1:end) = 4 / 3 * pdfb_thr(end-fssize+1:end);
    
    sorh        =    'h';                   % hard thresholding
    

    tvcx        =    wthresh(vcx, sorh, pdfb_thr);
    y           =    vec2pdfb(tvcx, s);
    %y           =    vec2pdfb(vcx, s);
  end
  
% $$$ function [y]=threshold_band(x,d,pfilt,dfilt,dlevs,rlevs)
% $$$ 
% $$$ for istr=1:length(dlevs)
% $$$   dstring(istr)=num2str(dlevs(istr));
% $$$ end
% $$$ 
% $$$ if ~exist(['nvar_ dstrin '.mat])
% $$$   nvar = pdfb_nest(size(d), pfilt, dfilt, dlevs);
% $$$   save(['nvar_ dstrin '.mat],nvar);
% $$$ else
% $$$    nvar =  load(['nvar_ dstrin '.mat],nvar);
% $$$ end
% $$$ 
% $$$   fac         =    abs(rlevs);
% $$$   
% $$$   pdfb_thr    =    fac * 3 *  nvar;
% $$$ 
% $$$   fssize      =    prod(s(end, :));	% finest scale size
% $$$ 
% $$$   pdfb_thr(end-fssize+1:end) = 4 / 3 * pdfb_thr(end-fssize+1:end);
% $$$ 
% $$$   sorh        =    'h';                   % hard thresholding
% $$$ 
% $$$   [vcz,s]     =    pdfb2vec(x);
% $$$   tvcz        =    wthresh(vcz, sorh, pdfb_thr);


  


