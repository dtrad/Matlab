function dyad_sig = paddyad(sig)
% paddyad -- Zero-fill signal to Dyadic length
%  Usage
%    dyad_sig = paddyad(sig)
%  Inputs
%    sig        a row or column vector
%  Outputs
%    dyad_sig   a vector of dyadic length
%               with contents taken from sig
%
%  See Also
%    cutdyad
%
	ss = size(sig);
	n = length(sig) ;
    J = ceil(log(n)/log(2));
	n1 = round(2^J);
    if n1 ~= n ,
      if ss(1) == n,
	  	dyad_sig = zeros(n1,1) ;
	  else
	    dyad_sig = zeros(1,n1) ;
	  end
	  dyad_sig(1:n) = sig(1:n);
    else
		dyad_sig = sig;
	end
    
%    
% Copyright (c) 1995, David L. Donoho
%
    
    
%   
% Part of WaveLab version .701
% Built Tuesday, January 30, 1996 8:25:59 PM
% This is Copyrighted Material
% For copying permissions see copying.m
% Comments? e-mail wavelab@playfair.stanford.edu
%   
    
   
% Auto name mapping to DOS conventions  Wednesday, January 31, 1996 5:36:19 PM
