function [xrec,wcx,twcx]=applyridget(x,perc,mask)
% apply the ridgelet transform, thresholding and inverse.
% usage  [xrec,wcx,twcx]=applyridget(x,perc,mask)
% mask is used to mute particular areas of the transform domain.
[nt nh]=size(x);

if (nargin<3) mask=ones(2*nt,2*nh);end
% Apply ridgelet transform

wcx=FastOrthoRidgeletTransform(x);

% do the thresholding 
%twcx =  HardThresh(wcx,T);

[T]=threshold(wcx,perc);
twcx=wcx;
twcx(abs(wcx)<T)=0;
%twcx=wcx;
% do the reconstruction
twcx=twcx.*mask;

xrec=real(Inv_FastOrthoRidgeletTrans(twcx));

return;



